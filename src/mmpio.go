// SPDX-License-Identifier: MIT
package main

import (
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	"sync"

	"gonum.org/v1/gonum/stat/distuv"
)

var outputPath string
var configPath string

type SumStatConf struct {
	Tag           string  `json:"tag"`
	Filepath      string  `json:"filepath"`
	ColChrom      string  `json:"col_chrom"`
	ColPos        string  `json:"col_pos"`
	ColRef        string  `json:"col_ref"`
	ColAlt        string  `json:"col_alt"`
	ColPval       string  `json:"col_pval"`
	ColBeta       string  `json:"col_beta"`
	ColSebeta     string  `json:"col_sebeta"`
	ColAf         string  `json:"col_af"`
	PvalThreshold float64 `json:"pval_threshold"`
	FmFilepath    string  `json:"fine_mapping_filepath"`
}

type HeterogeneityTest struct {
	Tag     string   `json:"tag"`
	Compare []string `json:"compare"`
}

type Conf struct {
	Inputs             []SumStatConf       `json:"inputs"`
	HeterogeneityTests []HeterogeneityTest `json:"heterogeneity_tests"`
}

type Cpra struct {
	Chrom string
	Pos   string
	Ref   string
	Alt   string
}

type Stats struct {
	Tag    string
	Pval   string
	Beta   string
	Sebeta string
	Af     string
}

type CpraStats struct {
	Cpra  Cpra
	Stats Stats
}

type FineMappingStats struct {
	Tag string
	Pip string
	Cs  string
}

type CpraFineMappingStats struct {
	Cpra             Cpra
	FineMappingStats FineMappingStats
}

type CombinedStats struct {
	Tag    string
	Pval   string
	Beta   string
	Sebeta string
	Af     string
	Pip    string
	Cs     string
}

func main() {
	conf := readConf(configPath)

	fmt.Println("[1/5] Checking variant selection...")
	selectedVariants := checkVariantSelection(conf.Inputs)

	fmt.Println("[2/5] Getting variant statistics...")
	statsVariants := findVariantStats(conf.Inputs, selectedVariants)

	fmt.Println("[3/5] Getting variant fine mapping statistics...")
	fmstatsVariants := findVariantFineMappingStats(conf.Inputs, selectedVariants)

	fmt.Println("[4/5] Combining variant GWAS statistics and fine mapping statistics...")
	combinedStatsVariants := combineStatsAndFmStats(statsVariants, fmstatsVariants)

	fmt.Printf("[5/5] Computing heterogeneity tests & writing output to %s ...\n", outputPath)
	writeMMPOutput(conf, combinedStatsVariants)
}

func init() {
	flag.StringVar(&outputPath, "output", "mmp.tsv", "Specify the output path (TSV)")
	flag.StringVar(&configPath, "config", "config.json", "Specify the configuration path (JSON)")
	flag.Parse()
}

func logCheck(message string, err error) {
	if err != nil {
		log.Fatal(":: ", message, " :: ", err)
	}
}

func readConf(filePath string) Conf {
	data, err := os.ReadFile(filePath)
	logCheck("reading configuration file", err)

	var conf Conf
	err = json.Unmarshal(data, &conf)
	logCheck("parsing JSON conf", err)

	// Validate JSON.
	// Go will not complain if there is a missing field in our input configuration file,
	// instead it will fill it with its default type value.
	// Since our fields are all required we manually check that all fields were provided
	// in the input configuration file.
	if conf.Inputs == nil {
		log.Fatal("Missing `inputs` field in the configuration file.")
	}
	if len(conf.Inputs) < 1 {
		log.Fatal("No summary stat provided in the configuration file. Need at least 1.")
	}
	for ii, input := range conf.Inputs {
		if input.Tag == "" {
			logMissingCol("tag", ii, "inputs")
		}
		if input.Filepath == "" {
			logMissingCol("filepath", ii, "inputs")
		}
		if input.ColChrom == "" {
			logMissingCol("col_chrom", ii, "inputs")
		}
		if input.ColPos == "" {
			logMissingCol("col_pos", ii, "inputs")
		}
		if input.ColRef == "" {
			logMissingCol("col_ref", ii, "inputs")
		}
		if input.ColAlt == "" {
			logMissingCol("col_alt", ii, "inputs")
		}
		if input.ColPval == "" {
			logMissingCol("col_pval", ii, "inputs")
		}
		if input.ColBeta == "" {
			logMissingCol("col_beta", ii, "inputs")
		}
		if input.ColSebeta == "" {
			logMissingCol("col_sebeta", ii, "inputs")
		}
		if input.ColAf == "" {
			logMissingCol("col_af", ii, "inputs")
		}
		if input.PvalThreshold == 0 {
			logMissingCol("pval_threshold", ii, "inputs")
		}
	}

	if conf.HeterogeneityTests == nil {
		log.Fatal("Missing `heterogeneity_tests` field in the configuration file.")
	}
	for jj, heterogeneity_test := range conf.HeterogeneityTests {
		if heterogeneity_test.Tag == "" {
			logMissingCol("tag", jj, "heterogeneity_tests")
		}
		if heterogeneity_test.Compare == nil {
			logMissingCol("compare", jj, "heterogeneity_tests")
		}
		if len(heterogeneity_test.Compare) < 2 {
			log.Fatal("Need at least 2 GWAS to run heterogeneity test. Instead got: ", heterogeneity_test.Compare)
		}
	}

	return conf
}

func logMissingCol(col_name string, element_index int, section string) {
	log.Fatal("Missing `", col_name, "` field of element #", element_index, " in the `", section, "` section of the configuration file. Check config.json.sample for reference.")
}

func checkVariantSelection(conf []SumStatConf) map[Cpra]bool {
	selectedVariants := make(map[Cpra]bool)

	var wg sync.WaitGroup
	ch := make(chan Cpra)

	for _, ssConf := range conf {
		wg.Add(1)
		go func(ssConf SumStatConf) {
			defer wg.Done()
			streamVariants(ssConf, ch)
		}(ssConf)
	}

	go func() {
		wg.Wait()
		close(ch)
	}()

	for cpra := range ch {
		selectedVariants[cpra] = true
	}

	return selectedVariants
}

func streamVariants(ssConf SumStatConf, ch chan<- Cpra) {
	fmt.Printf("- processing %s\n", ssConf.Tag)

	chRows := make(chan CpraStats)
	go readGzTsv(ssConf, chRows)

	for rec := range chRows {
		pvalUnparsed := rec.Stats.Pval
		pval, err := parseFloat64NaN(pvalUnparsed)
		logCheck("parsing p-value as float", err)

		if pval < ssConf.PvalThreshold {
			ch <- rec.Cpra
		}
	}

	fmt.Printf("* done %s\n", ssConf.Tag)
}

func findVariantStats(conf []SumStatConf, selectedVariants map[Cpra]bool) map[Cpra][]Stats {
	stats := make(map[Cpra][]Stats)

	var wg sync.WaitGroup
	ch := make(chan CpraStats)

	for _, ssConf := range conf {
		wg.Add(1)
		go func(ssConf SumStatConf) {
			defer wg.Done()
			streamStats(ssConf, selectedVariants, ch)
		}(ssConf)
	}

	go func() {
		wg.Wait()
		close(ch)
	}()

	for variant := range ch {
		multiStats, found := stats[variant.Cpra]
		if !found {
			var multiStats = []Stats{variant.Stats}
			stats[variant.Cpra] = multiStats
		} else {
			multiStats = append(multiStats, variant.Stats)
			stats[variant.Cpra] = multiStats
		}
	}

	return stats
}

func streamStats(ssConf SumStatConf, selectedVariants map[Cpra]bool, ch chan<- CpraStats) {
	fmt.Printf("- processing %s\n", ssConf.Tag)

	chRows := make(chan CpraStats)
	go readGzTsv(ssConf, chRows)

	for rec := range chRows {
		if _, found := selectedVariants[rec.Cpra]; found {
			ch <- rec
		}
	}

	fmt.Printf("* done %s\n", ssConf.Tag)
}

func parseCpra(record []string, fields map[string]int, ssConf SumStatConf) Cpra {
	fieldChrom, found := fields[ssConf.ColChrom]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColChrom, "` in header of input file `", ssConf.Filepath, "`.")
	}

	fieldPos, found := fields[ssConf.ColPos]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColPos, "` in header of input file `", ssConf.Filepath, "`.")
	}

	fieldRef, found := fields[ssConf.ColRef]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColRef, "` in header of input file `", ssConf.Filepath, "`.")
	}

	fieldAlt, found := fields[ssConf.ColAlt]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColAlt, "` in header of input file `", ssConf.Filepath, "`.")
	}

	return Cpra{
		record[fieldChrom],
		record[fieldPos],
		record[fieldRef],
		record[fieldAlt],
	}
}

func parseStats(record []string, fields map[string]int, ssConf SumStatConf) Stats {
	fieldPval, found := fields[ssConf.ColPval]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColPval, "` in header of input file `", ssConf.Filepath, "`.")
	}

	fieldBeta, found := fields[ssConf.ColBeta]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColBeta, "` in header of input file `", ssConf.Filepath, "`.")
	}

	fieldSebeta, found := fields[ssConf.ColSebeta]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColSebeta, "` in header of input file `", ssConf.Filepath, "`.")
	}

	fieldAf, found := fields[ssConf.ColAf]
	if !found {
		log.Fatal("Could not find column `", ssConf.ColAf, "` in header of input file `", ssConf.Filepath, "`.")
	}

	stats := Stats{
		ssConf.Tag,
		record[fieldPval],
		record[fieldBeta],
		record[fieldSebeta],
		record[fieldAf],
	}

	return stats
}

func sum(slice []float64) float64 {
	total := 0.0
	for _, v := range slice {
		total += v
	}
	return total
}

func writeMMPOutput(conf Conf, combinedStatsVariants map[Cpra][]CombinedStats) {
	var outRecords [][]string

	statsCols := []string{"pval", "beta", "sebeta", "af", "pip", "cs"}
	headerFields := []string{
		"chrom",
		"pos",
		"ref",
		"alt",
	}

	lenCpraFields := 4

	for _, ssConf := range conf.Inputs {
		for _, suffix := range statsCols {
			field := fmt.Sprintf("%s_%s", ssConf.Tag, suffix)
			headerFields = append(headerFields, field)
		}
	}

	// Loop to add meta fields for each heterogeneity test
	for _, test := range conf.HeterogeneityTests {
		headerFields = append(headerFields,
			fmt.Sprintf("%s_meta_beta", test.Tag),
			fmt.Sprintf("%s_meta_sebeta", test.Tag),
			fmt.Sprintf("%s_meta_pval", test.Tag),
			fmt.Sprintf("%s_meta_hetpval", test.Tag),
		)
	}

	outRecords = append(outRecords, headerFields)

	for cpra, cpraStats := range combinedStatsVariants {
		record := make([]string, len(headerFields))
		record[0] = cpra.Chrom
		record[1] = cpra.Pos
		record[2] = cpra.Ref
		record[3] = cpra.Alt

		for _, stats := range cpraStats {
			var offset int
			for ii, ssConf := range conf.Inputs {
				if ssConf.Tag == stats.Tag {
					offset = lenCpraFields + ii*len(statsCols)
				}
			}
			record[offset+0] = stats.Pval
			record[offset+1] = stats.Beta
			record[offset+2] = stats.Sebeta
			record[offset+3] = stats.Af
			record[offset+4] = stats.Pip
			record[offset+5] = stats.Cs
		}

		// Calculate meta stats here
		for _, test := range conf.HeterogeneityTests {
			var beta []float64
			var sebeta []float64
			for _, stats := range cpraStats {
				if contains(test.Compare, stats.Tag) {
					b, err := parseFloat64NaN(stats.Beta)
					logCheck("parsing beta as float", err)
					s, err := parseFloat64NaN(stats.Sebeta)
					logCheck("parsing sebeta as float", err)
					beta = append(beta, b)
					sebeta = append(sebeta, s)
				}
			}

			invVar := make([]float64, len(sebeta))
			for i := range invVar {
				invVar[i] = 1 / (sebeta[i] * sebeta[i])
			}
			effInvVar := make([]float64, len(beta))
			for i := range effInvVar {
				effInvVar[i] = beta[i] * invVar[i]
			}
			metaBeta := sum(effInvVar) / sum(invVar)

			metaSe := math.Sqrt(1 / sum(invVar))

			metaPVal := 2 * distuv.UnitNormal.Survival(math.Abs(sum(effInvVar))/math.Sqrt(sum(invVar)))

			// Calculate metaHetPVal here
			var betaDev []float64
			for i := range beta {
				betaDev = append(betaDev, invVar[i]*(beta[i]-metaBeta)*(beta[i]-metaBeta))
			}

			metaHetPVal := 1 - distuv.ChiSquared{
				K:   1,
				Src: nil,
			}.CDF(sum(betaDev))

			offset := lenCpraFields + len(conf.Inputs)*len(statsCols) + indexOfTest(test.Tag, conf.HeterogeneityTests)*len(statsCols)
			record[offset+0] = fmt.Sprintf("%f", metaBeta)
			record[offset+1] = fmt.Sprintf("%f", metaSe)
			record[offset+2] = fmt.Sprintf("%e", metaPVal)
			record[offset+3] = fmt.Sprintf("%e", metaHetPVal)
		}

		outRecords = append(outRecords, record)
	}

	outFile, err := os.Create(outputPath)
	logCheck("creating output file", err)
	defer outFile.Close()

	tsvWriter := csv.NewWriter(outFile)
	tsvWriter.Comma = '\t'
	tsvWriter.WriteAll(outRecords)
	err = tsvWriter.Error()
	logCheck("writing TSV output", err)
}

func contains(slice []string, item string) bool {
	for _, elem := range slice {
		if elem == item {
			return true
		}
	}
	return false
}

func indexOfTest(tag string, tests []HeterogeneityTest) int {
	for i, test := range tests {
		if test.Tag == tag {
			return i
		}
	}
	return -1
}

func readGzTsv(ssConf SumStatConf, ch chan<- CpraStats) {
	// Open gzip file for reading
	fReader, err := os.Open(ssConf.Filepath)
	logCheck("opening file", err)
	defer fReader.Close()

	gzReader, err := gzip.NewReader(fReader)
	logCheck("gunzip-ing file", err)
	defer gzReader.Close()

	// Parse as TSV
	tsvReader := csv.NewReader(gzReader)
	tsvReader.Comma = '\t'

	// Keep track of the TSV header
	rec, err := tsvReader.Read()
	logCheck("parsing TSV header", err)
	fields := make(map[string]int)
	for ii, field := range rec {
		fields[field] = ii
	}

	// Emit the TSV rows
	for {
		rec, err = tsvReader.Read()

		// Can't read data if end of file or parsing error
		if err == io.EOF {
			break
		}
		logCheck("parsing TSV row", err)

		// Emit the data if valid
		cpra := parseCpra(rec, fields, ssConf)
		stats := parseStats(rec, fields, ssConf)
		ch <- CpraStats{cpra, stats}
	}

	close(ch)
}

func parseFloat64NaN(input string) (float64, error) {
	const parseableNaN = "NaN"

	if input == "NA" {
		return strconv.ParseFloat(parseableNaN, 64)
	} else {
		return strconv.ParseFloat(input, 64)
	}
}

func findVariantFineMappingStats(conf []SumStatConf, selectedVariants map[Cpra]bool) map[Cpra][]FineMappingStats {
	fmstats := make(map[Cpra][]FineMappingStats)

	var wg sync.WaitGroup
	ch := make(chan CpraFineMappingStats)

	for _, ssConf := range conf {
		// Only add finemapping stats if a finemapping file was provided
		if ssConf.FmFilepath != "" {
			wg.Add(1)
			go func(ssConf SumStatConf) {
				defer wg.Done()
				streamFmStats(ssConf, selectedVariants, ch)
			}(ssConf)
		}
	}

	go func() {
		wg.Wait()
		close(ch)
	}()

	for variant := range ch {
		multiStats, found := fmstats[variant.Cpra]
		if !found {
			var multiStats = []FineMappingStats{variant.FineMappingStats}
			fmstats[variant.Cpra] = multiStats
		} else {
			multiStats = append(multiStats, variant.FineMappingStats)
			fmstats[variant.Cpra] = multiStats
		}
	}

	return fmstats
}

func streamFmStats(ssConf SumStatConf, selectedVariants map[Cpra]bool, ch chan<- CpraFineMappingStats) {
	fmt.Printf("- processing %s\n", ssConf.Tag)

	chRows := make(chan CpraFineMappingStats)
	go readTsv(ssConf, chRows)

	for rec := range chRows {
		if _, found := selectedVariants[rec.Cpra]; found {
			ch <- rec
		}
	}

	fmt.Printf("* done %s\n", ssConf.Tag)
}

func readTsv(ssConf SumStatConf, ch chan<- CpraFineMappingStats) {
	// Open file for reading
	fReader, err := os.Open(ssConf.FmFilepath)
	logCheck("opening file", err)
	defer fReader.Close()

	// Parse as TSV
	tsvReader := csv.NewReader(fReader)
	tsvReader.Comma = '\t'

	// Keep track of the TSV header
	rec, err := tsvReader.Read()
	logCheck("parsing TSV header", err)
	fields := make(map[string]int)
	for ii, field := range rec {
		fields[field] = ii
	}

	// Emit the TSV rows
	for {
		rec, err = tsvReader.Read()

		// Can't read data if end of file or parsing error
		if err == io.EOF {
			break
		}
		logCheck("parsing TSV row", err)

		// Emit the data if valid
		fmcpra := parseFmCpra(rec, fields)
		fmstats := parseFmStats(rec, fields, ssConf)
		ch <- CpraFineMappingStats{fmcpra, fmstats}
	}

	close(ch)
}

func parseFmCpra(record []string, fields map[string]int) Cpra {
	chromosomeWithPrefix := record[fields["chromosome"]]
	// Remove "chr" prefix
	chromosome := strings.TrimPrefix(chromosomeWithPrefix, "chr")
	position := record[fields["position"]]
	allele1 := record[fields["allele1"]]
	allele2 := record[fields["allele2"]]

	return Cpra{
		Chrom: chromosome,
		Pos:   position,
		Ref:   allele1,
		Alt:   allele2,
	}
}

func parseFmStats(record []string, fields map[string]int, ssConf SumStatConf) FineMappingStats {
	Pip := record[fields["cs_specific_prob"]]
	Cs := record[fields["cs"]]

	fmstats := FineMappingStats{
		ssConf.Tag,
		Pip,
		Cs,
	}

	return fmstats
}

func combineStatsAndFmStats(stats map[Cpra][]Stats, fmstats map[Cpra][]FineMappingStats) map[Cpra][]CombinedStats {
	combinedStatsVariants := make(map[Cpra][]CombinedStats)

	for cpra, statList := range stats {
		fmStatList := fmstats[cpra]

		combinedList := combineStatsAndFmStat(statList, fmStatList)
		combinedStatsVariants[cpra] = combinedList
	}

	return combinedStatsVariants
}

func combineStatsAndFmStat(statsList []Stats, fmStatList []FineMappingStats) []CombinedStats {
	combinedList := make([]CombinedStats, 0)

	for _, stat := range statsList {
		// Find matching FineMappingStats based on Tag
		var matchingFmStat FineMappingStats
		found := false
		for _, fmStat := range fmStatList {
			if fmStat.Tag == stat.Tag {
				matchingFmStat = fmStat
				found = true
				break
			}
		}

		// Combine Stats and FineMappingStats
		combined := CombinedStats{
			Tag:    stat.Tag,
			Pval:   stat.Pval,
			Beta:   stat.Beta,
			Sebeta: stat.Sebeta,
			Af:     stat.Af,
			Pip:    matchingFmStat.Pip,
			Cs:     matchingFmStat.Cs,
		}

		if !found {
			// Handle the case when no matching FineMappingStat is found
			combined.Pip = "NA" // or any other default value
			combined.Cs = "NA"  // or any other default value
		}

		combinedList = append(combinedList, combined)
	}

	return combinedList
}
