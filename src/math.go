package main

import (
	"math"

	"gonum.org/v1/gonum/stat/distuv"
)

// For heterogeneity test
type OutputMetaStats struct {
	Beta    string
	SEBeta  string
	PVal    string
	HetPVal string
}

func ComputeHeterogeneityTest(Betas []float64, SEBetas []float64) OutputMetaStats {
	invVar := make([]float64, len(SEBetas))
	for i := range invVar {
		invVar[i] = 1 / (SEBetas[i] * SEBetas[i])
	}

	effInvVar := make([]float64, len(Betas))
	for i := range effInvVar {
		effInvVar[i] = Betas[i] * invVar[i]
	}

	metaBeta := sum(effInvVar) / sum(invVar)
	metaSEBeta := math.Sqrt(1 / sum(invVar))
	metaPVal := 2 * distuv.UnitNormal.Survival(math.Abs(sum(effInvVar))/math.Sqrt(sum(invVar)))

	// Calculate metaHetPVal here
	var betaDev []float64
	for i := range Betas {
		betaDev = append(betaDev, invVar[i]*(Betas[i]-metaBeta)*(Betas[i]-metaBeta))
	}

	metaHetPVal := 1 - distuv.ChiSquared{
		K:   1,
		Src: nil,
	}.CDF(sum(betaDev))

	// Convert values to string for outputting and return
	return OutputMetaStats{
		Beta:    formatFloat(metaBeta),
		SEBeta:  formatFloat(metaSEBeta),
		PVal:    formatFloat(metaPVal),
		HetPVal: formatFloat(metaHetPVal),
	}
}
