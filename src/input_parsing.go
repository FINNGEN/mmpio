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
)

type CPRA struct {
	Chrom string
	Pos   string
	Ref   string
	Alt   string
}

type SummaryStats struct {
	PVal   string
	Beta   string
	SEBeta string
	AF     string
}

// This is using struct embedding, see https://gobyexample.com/struct-embedding
// It is useful like this because in some places we need the CPRA struct without
// the stats attached.
type InputSummaryStatsRow struct {
	Tag string
	CPRA
	SummaryStats
}

type InputFinemapRow struct {
	Tag string
	CPRA
	PIP string
	CS  string
}

type OutputStats struct {
	Tag    string
	PVal   string
	Beta   string
	SEBeta string
	AF     string
	PIP    string
	CS     string
}

func streamVariantsAboveThreshold(inputConf InputConf, cpraChannel chan<- CPRA) {
	fmt.Printf("- processing %s\n", inputConf.Tag)

	parsedRowChannel := make(chan InputSummaryStatsRow)
	go streamSummaryStatsFile(inputConf, parsedRowChannel)

	for row := range parsedRowChannel {
		parsedPVal, err := parseFloat64NaN(row.PVal)
		logCheck("parsing p-value as float", err)

		if parsedPVal < inputConf.PValThreshold {
			cpraChannel <- row.CPRA
		}
	}

	fmt.Printf("* done %s\n", inputConf.Tag)
}

func streamRowsFromSelection(inputConf InputConf, selectedVariants map[CPRA]bool, selectedRowChannel chan<- InputSummaryStatsRow) {
	fmt.Printf("- processing %s\n", inputConf.Tag)

	parsedRowChannel := make(chan InputSummaryStatsRow)
	go streamSummaryStatsFile(inputConf, parsedRowChannel)

	for row := range parsedRowChannel {
		if _, found := selectedVariants[row.CPRA]; found {
			selectedRowChannel <- row
		}
	}

	fmt.Printf("* done %s\n", inputConf.Tag)
}

// TODO(future)::STREAM-STRUCT is there a way to pass a mapping of columns (string) to struct,
// so that we stream any tabular file to a chan of struct (struct being any).
// maybe with generics?
// or with struct tags? https://go.dev/wiki/Well-known-struct-tags
func streamSummaryStatsFile(inputConf InputConf, parsedRowChannel chan<- InputSummaryStatsRow) {
	rowChannel := make(chan map[string]string)
	go streamTsv(inputConf.Filepath, "gzip", rowChannel)

	for row := range rowChannel {
		chrom, found := row[inputConf.ColChrom]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColChrom, "` in header of input file `", inputConf.Filepath, "`.")
		}

		pos, found := row[inputConf.ColPos]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColPos, "` in header of input file `", inputConf.Filepath, "`.")
		}

		ref, found := row[inputConf.ColRef]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColRef, "` in header of input file `", inputConf.Filepath, "`.")
		}

		alt, found := row[inputConf.ColAlt]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColAlt, "` in header of input file `", inputConf.Filepath, "`.")
		}

		pval, found := row[inputConf.ColPVal]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColPVal, "` in header of input file `", inputConf.Filepath, "`.")
		}

		beta, found := row[inputConf.ColBeta]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColBeta, "` in header of input file `", inputConf.Filepath, "`.")
		}

		seBeta, found := row[inputConf.ColSEBeta]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColSEBeta, "` in header of input file `", inputConf.Filepath, "`.")
		}

		af, found := row[inputConf.ColAF]
		if !found {
			log.Fatal("Could not find column `", inputConf.ColAF, "` in header of input file `", inputConf.Filepath, "`.")
		}

		parsedRow := InputSummaryStatsRow{
			Tag:          inputConf.Tag,
			CPRA:         CPRA{chrom, pos, ref, alt},
			SummaryStats: SummaryStats{pval, beta, seBeta, af},
		}

		parsedRowChannel <- parsedRow
	}
	close(parsedRowChannel)
}

// TODO(future) see ::STREAM-STRUCT
func streamFinemapFile(inputConf InputConf, parsedRowChannel chan<- InputFinemapRow) {
	colCPRA := "v"
	colPIP := "cs_specific_prob"
	colCS := "cs"

	fmt.Printf("- processing %s\n", inputConf.Tag)

	rowChannel := make(chan map[string]string)
	go streamTsv(inputConf.FinemapFilepath, "uncompressed", rowChannel)

	for row := range rowChannel {
		// Parse the CPRA from assumed "C:P:R:A" format
		cpra, found := row[colCPRA]
		if !found {
			log.Fatal("Could not find column `", colCPRA, "` in header of input file `", inputConf.FinemapFilepath, "`.")
		}
		splitCPRA := strings.Split(cpra, ":")
		if len(splitCPRA) != 4 {
			log.Fatal("Could not parse CPRA from value `", cpra, "`.")
		}
		chrom := splitCPRA[0]
		pos := splitCPRA[1]
		ref := splitCPRA[2]
		alt := splitCPRA[3]

		pip, found := row[colPIP]
		if !found {
			log.Fatal("Could not find column `", colPIP, "` in header of input file `", inputConf.FinemapFilepath, "`.")
		}

		cs, found := row[colCS]
		if !found {
			log.Fatal("Could not find column `", colCS, "` in header of input file `", inputConf.FinemapFilepath, "`.")
		}

		parsedRow := InputFinemapRow{
			Tag:  inputConf.Tag,
			CPRA: CPRA{chrom, pos, ref, alt},
			PIP:  pip,
			CS:   cs,
		}

		parsedRowChannel <- parsedRow
	}

	fmt.Printf("* done %s\n", inputConf.Tag)
}

func streamTsv(filepath string, compressionType string, rowChannel chan<- map[string]string) {
	// For now just assume we only deal with TSV files.
	// Don't implement support for other input formats until we need it.
	//
	// Given an input file like this:
	//     col1  col2  col3
	//       a1    a2    a3
	//       b1    b2    b3
	// Send this data over the channel:
	//     map[string]string{"col1": "a1", "col2": "a2", "col3": "a3"}
	//     map[string]string{"col1": "b1", "col2": "b2", "col3": "b3"}
	//     map[string]string{"col1": "c1", "col2": "c2", "col3": "c3"}

	// Open file for reading
	fReader, err := os.Open(filepath)
	logCheck("opening file", err)
	defer fReader.Close()

	// Uncompress the file if necessary
	var dataReader io.Reader

	switch compressionType {
	case "uncompressed":
		dataReader = fReader

	case "gzip":
		gzReader, err := gzip.NewReader(fReader)
		logCheck("gunzip-ing file", err)
		defer gzReader.Close()
		dataReader = gzReader

	default:
		log.Fatal("Unrecognized compression type `", compressionType, "`. Possible values are: uncompressed, gzip.")
	}

	// Parse as TSV
	tsvReader := csv.NewReader(dataReader)
	tsvReader.Comma = '\t'

	// Keep track of the TSV header
	header, err := tsvReader.Read()
	logCheck("parsing TSV header", err)

	// Emit the rows over the channel
	for {
		row, err := tsvReader.Read()

		// Can't read more data if end of file or parsing error
		if err == io.EOF {
			break
		}
		logCheck("parsing TSV row", err)

		// If valid, emit the data as a map of Header -> Row value
		rowWithHeader := make(map[string]string)
		for ii := 0; ii < len(header); ii++ {
			rowWithHeader[header[ii]] = row[ii]
		}
		rowChannel <- rowWithHeader
	}

	close(rowChannel)
}
