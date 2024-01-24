// SPDX-License-Identifier: MIT

// TODO
// - Refactor writing to a TSV, we don't need to deal with column indices.
//   See related ::STREAM-STRUCT

package main

import (
	"encoding/csv"
	"fmt"
	"os"
)

func writeMMPOutput(conf Conf, combinedStatsVariants map[CPRA][]OutputStats) {
	var outRecords [][]string

	statsCols := []string{"pval", "beta", "sebeta", "af", "pip", "cs"}
	headerFields := []string{
		"chrom",
		"pos",
		"ref",
		"alt",
	}

	lenCpraFields := 4

	for _, inputConf := range conf.Inputs {
		for _, suffix := range statsCols {
			field := fmt.Sprintf("%s_%s", inputConf.Tag, suffix)
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

	for cpra, multipleStats := range combinedStatsVariants {
		// Initialize the record
		record := make([]string, len(headerFields))
		record[0] = cpra.Chrom
		record[1] = cpra.Pos
		record[2] = cpra.Ref
		record[3] = cpra.Alt

		for ii := 4; ii < len(headerFields); ii++ {
			// If a summary stats file doesn't contain a given CPRA, then
			// we will show "NA" in the output for its stats.
			// If has the given CPRA, then the "NA" will be overwritten
			// by the actual summary stats values in the next step.
			record[ii] = outputDefaultMissingValue
		}

		// Add summary statistics for each of the input
		for _, stats := range multipleStats {
			var offset int
			for ii, inputConf := range conf.Inputs {
				if inputConf.Tag == stats.Tag {
					offset = lenCpraFields + ii*len(statsCols)
				}
			}
			record[offset+0] = stats.PVal
			record[offset+1] = stats.Beta
			record[offset+2] = stats.SEBeta
			record[offset+3] = stats.AF
			record[offset+4] = stats.PIP
			record[offset+5] = stats.CS
		}

		// Check tags with stats for het test
		tagsWithStats := make(map[string]bool)
		for _, stats := range multipleStats {
			if stats.Beta != "NA" && stats.SEBeta != "NA" {
				tagsWithStats[stats.Tag] = true
			}
		}

		// Calculate meta stats here
		for _, test := range conf.HeterogeneityTests {
			// Check the test has necessary data
			hasNecessaryData := true
			for _, tagCompare := range test.Compare {
				_, found := tagsWithStats[tagCompare]
				if !found {
					hasNecessaryData = false
					break
				}
			}

			var metaStats OutputMetaStats
			if hasNecessaryData {
				var betas []float64
				var sebetas []float64
				for _, stats := range multipleStats {
					if contains(test.Compare, stats.Tag) {
						beta, err := parseFloat64NaN(stats.Beta)
						logCheck("parsing beta as float", err)
						betas = append(betas, beta)

						sebeta, err := parseFloat64NaN(stats.SEBeta)
						logCheck("parsing sebeta as float", err)
						sebetas = append(sebetas, sebeta)
					}
				}
				metaStats = ComputeHeterogeneityTest(betas, sebetas)
			} else {
				// Don't compute the meta stats if some stats are missing
				metaStats = OutputMetaStats{
					Beta:    "NA",
					SEBeta:  "NA",
					PVal:    "NA",
					HetPVal: "NA",
				}
			}

			offset := lenCpraFields + len(conf.Inputs)*len(statsCols) + indexOfTest(test.Tag, conf.HeterogeneityTests)*len(statsCols)
			record[offset+0] = metaStats.Beta
			record[offset+1] = metaStats.SEBeta
			record[offset+2] = metaStats.PVal
			record[offset+3] = metaStats.HetPVal
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

func indexOfTest(tag string, tests []HeterogeneityTestConf) int {
	for i, test := range tests {
		if test.Tag == tag {
			return i
		}
	}
	return -1
}
