// SPDX-License-Identifier: MIT
package main

import (
	"fmt"
)

func main() {
	conf := readConf(configPath)

	fmt.Println("[1/5] Scanning input files for variant selection...")
	selectedVariants := scanForVariantSelection(conf.Inputs)

	fmt.Println("[2/5] Finding variant statistics based on the variant selection...")
	statsVariants := findVariantStats(conf.Inputs, selectedVariants)

	fmt.Println("[3/5] Finding variant fine mapping statistics...")
	fmstatsVariants := findVariantFinemapping(conf.Inputs, selectedVariants)

	// TODO remove with refactoring to same struct
	fmt.Println("[4/5] Combining variant GWAS statistics and fine mapping statistics...")
	combinedStatsVariants := combineStatsAndFmStats(statsVariants, fmstatsVariants)

	fmt.Printf("[5/5] Computing heterogeneity tests & writing output to %s ...\n", outputPath)
	writeMMPOutput(conf, combinedStatsVariants)
}
