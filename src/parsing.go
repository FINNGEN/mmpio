// SPDX-License-Identifier: MIT
package main

import (
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"sync"
)

// Structs related to the input files
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

// TODO this should just be changed to have the finemapping fields in the same struct as pval, beta, etc.
// We want:
// - input tag:A :  CPRA => pval, beta , se beta, af, pip, cs
// - input tag:B :  CPRA => pval, beta , se beta, af, pip, cs
// ...
type CpraFineMappingStats struct {
	Cpra             Cpra
	FineMappingStats FineMappingStats
}

// TODO either use this and skip 'Stats' and 'CpraFineMappingStats', or remove it in favor of extending 'Stats' with finemapping fields
type CombinedStats struct {
	Tag    string
	Pval   string
	Beta   string
	Sebeta string
	Af     string
	Pip    string
	Cs     string
}

func scanForVariantSelection(conf []SumStatConf) map[Cpra]bool {
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

func findVariantFinemapping(conf []SumStatConf, selectedVariants map[Cpra]bool) map[Cpra][]FineMappingStats {
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
