// SPDX-License-Identifier: MIT
package main

import (
	"fmt"
	"sync"
)

const outputDefaultMissingValue = "NA"

func main() {
	cliInit()
	conf := readConf(configPath)

	fmt.Println("[1/4] Scanning input files for variant selection...")
	selectedVariants := scanForVariantSelection(conf)

	fmt.Println("[2/4] Finding variant statistics based on the variant selection...")
	variantStats := findVariantStats(conf, selectedVariants)

	fmt.Println("[3/4] Combining finemapping statistics...")
	combineFinemapping(conf, variantStats)

	fmt.Printf("[4/4] Computing heterogeneity tests & writing output to %s ...\n", outputPath)
	writeMMPOutput(conf, variantStats)
}

func scanForVariantSelection(conf Conf) map[CPRA]bool {
	selectedVariants := make(map[CPRA]bool)

	var wg sync.WaitGroup
	cpraChannel := make(chan CPRA)

	for _, inputConf := range conf.Inputs {
		wg.Add(1)
		go func(inputConf InputConf) {
			defer wg.Done()
			streamVariantsAboveThreshold(inputConf, cpraChannel)
		}(inputConf)
	}

	go func() {
		wg.Wait()
		close(cpraChannel)
	}()

	for cpra := range cpraChannel {
		selectedVariants[cpra] = true
	}

	return selectedVariants
}

func findVariantStats(conf Conf, selectedVariants map[CPRA]bool) map[CPRA][]OutputStats {
	variantMultipleStats := make(map[CPRA][]OutputStats)

	var wg sync.WaitGroup
	selectedRowChannel := make(chan InputSummaryStatsRow)

	for _, inputConf := range conf.Inputs {
		wg.Add(1)
		go func(inputConf InputConf) {
			defer wg.Done()
			streamRowsFromSelection(inputConf, selectedVariants, selectedRowChannel)
		}(inputConf)
	}

	go func() {
		wg.Wait()
		close(selectedRowChannel)
	}()

	for parsedRow := range selectedRowChannel {
		outputStats := OutputStats{
			Tag:    parsedRow.Tag,
			PVal:   parsedRow.PVal,
			Beta:   parsedRow.Beta,
			SEBeta: parsedRow.SEBeta,
			AF:     parsedRow.AF,

			// These will be eventually filled with the finemapping values,
			// if a finemapping file was provided for this input.
			PIP: outputDefaultMissingValue,
			CS:  outputDefaultMissingValue,
		}

		multipleOutputStats, found := variantMultipleStats[parsedRow.CPRA]
		if !found {
			var multipleOutputStats = []OutputStats{outputStats}
			variantMultipleStats[parsedRow.CPRA] = multipleOutputStats
		} else {
			multipleOutputStats = append(multipleOutputStats, outputStats)
			variantMultipleStats[parsedRow.CPRA] = multipleOutputStats
		}
	}

	return variantMultipleStats
}

func combineFinemapping(conf Conf, variantStats map[CPRA][]OutputStats) {
	// We need this:
	// Tag => CPRA => InputFinemapRow
	// that is gathered by reading the finemap files
	//
	// Then we iterate over variantStats,
	// for each CPRA,
	// we look at each Tag,
	// and if there is  Tag => CPRA  match in the above, then we add
	// the finemap stats.
	finemapStatsGathering := make(map[string]map[CPRA]InputFinemapRow)

	var wg sync.WaitGroup
	finemapRowChannel := make(chan InputFinemapRow)

	for _, inputConf := range conf.Inputs {
		if inputConf.FinemapFilepath != "" {
			wg.Add(1)
			go func(inputConf InputConf) {
				defer wg.Done()
				streamFinemapFile(inputConf, finemapRowChannel)
			}(inputConf)
		}
	}

	go func() {
		wg.Wait()
		close(finemapRowChannel)
	}()

	for finemapRow := range finemapRowChannel {
		tagData, found := finemapStatsGathering[finemapRow.Tag]
		if !found {
			tagData = make(map[CPRA]InputFinemapRow)
		}

		tagData[finemapRow.CPRA] = finemapRow
		finemapStatsGathering[finemapRow.Tag] = tagData
	}

	// Now time to combine with variantStats
	for cpra, multipleOutputStats := range variantStats {
		for idxTag, outputStats := range multipleOutputStats {
			if _, tagFound := finemapStatsGathering[outputStats.Tag]; tagFound {
				if finemapStats, cpraFound := finemapStatsGathering[outputStats.Tag][cpra]; cpraFound {
					variantStats[cpra][idxTag].PIP = finemapStats.PIP
					variantStats[cpra][idxTag].CS = finemapStats.CS
				}
			}
		}
	}
}
