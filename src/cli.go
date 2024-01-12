// SPDX-License-Identifier: MIT
package main

import (
	"encoding/json"
	"flag"
	"log"
	"os"
)

var outputPath string
var configPath string

// Structs relating to configuration.
// Basically a 1:1 mapping to the conf fields and tree structure.
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
	FmFilepath    string  `json:"fine_mapping_filepath"` // TODO can we make this true optional, as in: allow this field to be omitted completely from the config file
}

type HeterogeneityTest struct {
	Tag     string   `json:"tag"`
	Compare []string `json:"compare"`
}

type Conf struct {
	Inputs             []SumStatConf       `json:"inputs"`
	HeterogeneityTests []HeterogeneityTest `json:"heterogeneity_tests"`
}

func init() {
	flag.StringVar(&outputPath, "output", "mmp.tsv", "Specify the output path (TSV)")
	flag.StringVar(&configPath, "config", "config.json", "Specify the configuration path (JSON)")
	flag.Parse()
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
	// TODO make sure we allow for the fimemapping_path conf field to not exist.
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
