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
	rowChannel := make(chan []string)
	requestedColumns := []string{
		inputConf.ColChrom,
		inputConf.ColPos,
		inputConf.ColRef,
		inputConf.ColAlt,
		inputConf.ColPVal,
		inputConf.ColBeta,
		inputConf.ColSEBeta,
		inputConf.ColAF,
	}
	go streamTsv(inputConf.Filepath, "gzip", requestedColumns, rowChannel)

	for row := range rowChannel {
		chrom := row[0]
		pos := row[1]
		ref := row[2]
		alt := row[3]
		pval := row[4]
		beta := row[5]
		seBeta := row[6]
		af := row[7]

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

	rowChannel := make(chan []string)
	requestedColumns := []string{
		colCPRA,
		colPIP,
		colCS,
	}
	go streamTsv(inputConf.FinemapFilepath, "uncompressed", requestedColumns, rowChannel)

	for row := range rowChannel {
		cpra := row[0]
		pip := row[1]
		cs := row[2]

		// Parse the CPRA from assumed "C:P:R:A" format
		splitCPRA := strings.Split(cpra, ":")
		if len(splitCPRA) != 4 {
			log.Fatal("Could not parse CPRA from value `", cpra, "`.")
		}
		chrom := splitCPRA[0]
		pos := splitCPRA[1]
		ref := splitCPRA[2]
		alt := splitCPRA[3]

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

func streamTsv(filepath string, compressionType string, columns []string, rowChannel chan<- []string) {
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

	headerToIndex := make(map[string]int)
	for ii, headerColumn := range header {
		headerToIndex[headerColumn] = ii
	}

	// Derive the field indices we want from the header
	requestedColIndices := make([]int, len(columns))
	for ii, requestedColumn := range columns {
		headerColumnIndex, found := headerToIndex[requestedColumn]
		if found {
			requestedColIndices[ii] = headerColumnIndex
		} else {
			log.Fatal("Could not find column `", requestedColumn, "` in header of input file `", filepath, "`. Header: ", header)
		}
	}

	// Emit the rows over the channel
	for {
		row, err := tsvReader.Read()

		// Can't read more data if end of file or parsing error
		if err == io.EOF {
			break
		}
		logCheck("parsing TSV row", err)

		// This variable needs to be initialized *inside* the for loop.
		// If it is assigned outside (in the hope of getting some performance improvements?),
		// then the slice will be concurrently read and written, causing bad data parsing down the line.
		// This was also caught by the go data race detector.
		rowFromColumns := make([]string, len(columns))
		for ii, requestedColIndex := range requestedColIndices {
			rowFromColumns[ii] = row[requestedColIndex]
		}

		rowChannel <- rowFromColumns
	}

	close(rowChannel)
}
