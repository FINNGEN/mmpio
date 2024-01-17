// SPDX-License-Identifier: MIT

// TODO
// FIN: - bug (kind of) in that the data output has fields with "" instead of "NA" when a variant
// is present in one dataset but not in the other.
// - Also, het test values should be NA when only one dataset has values (currently has values and p=1)

// Refactor writing to a TSV, we don't need to deal with column indices

// Refactor the het test out of the TSV write?

package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"

	"gonum.org/v1/gonum/stat/distuv"
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
			record[offset+0] = stats.PVal
			record[offset+1] = stats.Beta
			record[offset+2] = stats.SEBeta
			record[offset+3] = stats.AF
			record[offset+4] = stats.PIP
			record[offset+5] = stats.CS
		}

		// Calculate meta stats here
		for _, test := range conf.HeterogeneityTests {
			var beta []float64
			var sebeta []float64
			for _, stats := range cpraStats {
				if contains(test.Compare, stats.Tag) {
					b, err := parseFloat64NaN(stats.Beta)
					logCheck("parsing beta as float", err)
					s, err := parseFloat64NaN(stats.SEBeta)
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

func indexOfTest(tag string, tests []HeterogeneityTestConf) int {
	for i, test := range tests {
		if test.Tag == tag {
			return i
		}
	}
	return -1
}
