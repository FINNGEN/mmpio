// SPDX-License-Identifier: MIT
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
)

var outputPath string
var configPath string
var showVersion bool

// Get the program version from git.
// This should be passed as a build time variable, for example:
// go build -ldflags "-X main.MMPioVersion=$(git describe --tags)"
var MMPioVersion string

type InputConf struct {
	Tag             string  `json:"tag"`
	Filepath        string  `json:"filepath"`
	ColChrom        string  `json:"col_chrom"`
	ColPos          string  `json:"col_pos"`
	ColRef          string  `json:"col_ref"`
	ColAlt          string  `json:"col_alt"`
	ColPVal         string  `json:"col_pval"`
	ColBeta         string  `json:"col_beta"`
	ColSEBeta       string  `json:"col_sebeta"`
	ColAF           string  `json:"col_af"`
	PValThreshold   float64 `json:"pval_threshold"`
	FinemapFilepath string  `json:"finemap_filepath"`
}

type HeterogeneityTestConf struct {
	Tag     string   `json:"tag"`
	Compare []string `json:"compare"`
}

type Conf struct {
	Inputs             []InputConf             `json:"inputs"`
	HeterogeneityTests []HeterogeneityTestConf `json:"heterogeneity_tests"`
}

func init() {
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage of %s (%s):\n", os.Args[0], MMPioVersion)
		flag.PrintDefaults()
	}
	flag.StringVar(&configPath, "config", "config.json", "Specify the configuration path (JSON)")
	flag.StringVar(&outputPath, "output", "mmp.tsv", "Specify the output path (TSV)")

	flag.BoolVar(&showVersion, "version", false, "Show MMP::io version")
	flag.Parse()

	if showVersion {
		fmt.Fprintf(flag.CommandLine.Output(), "%s\n", MMPioVersion)
		os.Exit(0)
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
			logMissingKey("tag", ii, "inputs")
		}
		if input.Filepath == "" {
			logMissingKey("filepath", ii, "inputs")
		}
		if input.ColChrom == "" {
			logMissingKey("col_chrom", ii, "inputs")
		}
		if input.ColPos == "" {
			logMissingKey("col_pos", ii, "inputs")
		}
		if input.ColRef == "" {
			logMissingKey("col_ref", ii, "inputs")
		}
		if input.ColAlt == "" {
			logMissingKey("col_alt", ii, "inputs")
		}
		if input.ColPVal == "" {
			logMissingKey("col_pval", ii, "inputs")
		}
		if input.ColBeta == "" {
			logMissingKey("col_beta", ii, "inputs")
		}
		if input.ColSEBeta == "" {
			logMissingKey("col_sebeta", ii, "inputs")
		}
		if input.ColAF == "" {
			logMissingKey("col_af", ii, "inputs")
		}
		if input.PValThreshold == 0 {
			logMissingKey("pval_threshold", ii, "inputs")
		}
		// We don't check for the "fine_mapping_path" configuration key as it is optional.
	}

	if conf.HeterogeneityTests == nil {
		log.Fatal("Missing `heterogeneity_tests` field in the configuration file.")
	}
	for jj, heterogeneity_test := range conf.HeterogeneityTests {
		if heterogeneity_test.Tag == "" {
			logMissingKey("tag", jj, "heterogeneity_tests")
		}
		if heterogeneity_test.Compare == nil {
			logMissingKey("compare", jj, "heterogeneity_tests")
		}
		if len(heterogeneity_test.Compare) < 2 {
			log.Fatal("Need at least 2 GWAS to run heterogeneity test. Instead got: ", heterogeneity_test.Compare)
		}
	}

	return conf
}

func logMissingKey(col_name string, element_index int, section string) {
	log.Fatal("Missing `", col_name, "` key of element #", element_index, " in the `", section, "` section of the configuration file. Check config.json.sample for reference.")
}
