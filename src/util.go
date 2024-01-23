// SPDX-License-Identifier: MIT
package main

import (
	"log"
	"strconv"
)

func logCheck(message string, err error) {
	if err != nil {
		log.Fatal(":: ", message, " :: ", err)
	}
}

func sum(slice []float64) float64 {
	total := 0.0
	for _, v := range slice {
		total += v
	}
	return total
}

func contains(slice []string, item string) bool {
	for _, elem := range slice {
		if elem == item {
			return true
		}
	}
	return false
}

func containsNaN(slice []float64) bool {
	floatNaN, _ := strconv.ParseFloat("NaN", 64)

	for _, elem := range slice {
		if elem == floatNaN {
			return true
		}
	}

	return false
}

func parseFloat64NaN(input string) (float64, error) {
	const parseableNaN = "NaN"

	if input == "NA" {
		return strconv.ParseFloat(parseableNaN, 64)
	} else {
		return strconv.ParseFloat(input, 64)
	}
}

func formatFloat(number float64) string {
	var withDecimalExponent byte = 'e'
	precisionExactSmallest := -1

	return strconv.FormatFloat(number, withDecimalExponent, precisionExactSmallest, 64)
}
