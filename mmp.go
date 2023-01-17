// SPDX-License-Identifier: MIT
package main

import (
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"sync"
)

const OutputPath = "mmp.tsv"

type SumStatConf struct {
	Tag           string  `json:"tag"`
	Filepath      string  `json:"filepath"`
	ColChrom      string  `json:"col_chrom"`
	ColPos        string  `json:"col_pos"`
	ColRef        string  `json:"col_ref"`
	ColAlt        string  `json:"col_alt"`
	ColPval       string  `json:"col_pval"`
	ColBeta       string  `json:"col_beta"`
	ColAf         string  `json:"col_af"`
	PvalThreshold float64 `json:"pval_threshold"`
}

type Cpra struct {
	Chrom string
	Pos   string
	Ref   string
	Alt   string
}

type Stats struct {
	Tag  string
	Pval string
	Beta string
	Af   string
}

type CpraStats struct {
	Cpra  Cpra
	Stats Stats
}

func main() {
	conf := readConf("config.json")

	fmt.Println("[1/3] Checking variant selection...")
	selectedVariants := checkVariantSelection(conf)

	fmt.Println("[2/3] Getting variant statistics...")
	statsVariants := findVariantStats(conf, selectedVariants)

	fmt.Printf("[3/3] Writing output to %s ...\n", OutputPath)
	writeMMPOutput(conf, statsVariants)
}

func logCheck(message string, err error) {
	if err != nil {
		log.Fatal(":: ", message, " :: ", err)
	}
}

func readConf(filePath string) []SumStatConf {
	data, err := os.ReadFile(filePath)
	logCheck("reading configuration file", err)

	var conf []SumStatConf
	err = json.Unmarshal(data, &conf)
	logCheck("parsing JSON conf", err)

	return conf
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
		pval, err := strconv.ParseFloat(pvalUnparsed, 64)
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
	return Cpra{
		record[fields[ssConf.ColChrom]],
		record[fields[ssConf.ColPos]],
		record[fields[ssConf.ColRef]],
		record[fields[ssConf.ColAlt]],
	}
}

func parseStats(record []string, fields map[string]int, ssConf SumStatConf) Stats {
	pval := record[fields[ssConf.ColPval]]
	beta := record[fields[ssConf.ColBeta]]
	af := record[fields[ssConf.ColAf]]

	stats := Stats{
		ssConf.Tag,
		pval,
		beta,
		af,
	}

	return stats
}

func writeMMPOutput(conf []SumStatConf, statsVariants map[Cpra][]Stats) {
	var outRecords [][]string

	statsCols := []string{"pval", "beta", "af"}
	headerFields := []string{
		"chrom",
		"pos",
		"ref",
		"alt",
	}
	lenCpraFields := 4

	for _, ssConf := range conf {
		for _, suffix := range statsCols {
			field := fmt.Sprintf("%s_%s", ssConf.Tag, suffix)
			headerFields = append(headerFields, field)
		}
	}
	outRecords = append(outRecords, headerFields)

	for cpra, cpraStats := range statsVariants {
		record := make([]string, len(headerFields))
		record[0] = cpra.Chrom
		record[1] = cpra.Pos
		record[2] = cpra.Ref
		record[3] = cpra.Alt

		for _, stats := range cpraStats {
			var offset int
			for ii, ssConf := range conf {
				if ssConf.Tag == stats.Tag {
					offset = lenCpraFields + ii*len(statsCols)
				}
			}
			record[offset+0] = stats.Pval
			record[offset+1] = stats.Beta
			record[offset+2] = stats.Af
		}
		outRecords = append(outRecords, record)
	}

	outFile, err := os.Create(OutputPath)
	logCheck("creating output file", err)
	defer outFile.Close()

	tsvWriter := csv.NewWriter(outFile)
	tsvWriter.Comma = '\t'
	tsvWriter.WriteAll(outRecords)
	err = tsvWriter.Error()
	logCheck("writing TSV output", err)
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
