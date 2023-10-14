// R_init_DECIPHER.c

void R_init_DECIPHER(DllInfo *info);

// ConsensusSequence.c

SEXP consensusSequence(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps);

SEXP consensusSequenceAA(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps);

SEXP consensusProfile(SEXP x, SEXP weight, SEXP structs);

SEXP consensusProfileAA(SEXP x, SEXP weight, SEXP structs);

SEXP colScores(SEXP x, SEXP subset, SEXP subMatrix, SEXP go, SEXP ge, SEXP terminalGaps, SEXP weights, SEXP structs, SEXP dbnMatrix);

SEXP colScoresAA(SEXP x, SEXP subset, SEXP subMatrix, SEXP go, SEXP ge, SEXP terminalGaps, SEXP weights, SEXP structs, SEXP hecMatrix);

SEXP shiftGaps(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP gl, SEXP sc, SEXP thresh, SEXP weights);

SEXP shiftGapsAA(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP gl, SEXP sc, SEXP thresh, SEXP weights);

// DistanceMatrix.c

SEXP distMatrix(SEXP x, SEXP t, SEXP terminalGaps, SEXP penalizeGapLetters, SEXP fullMatrix, SEXP output, SEXP e, SEXP minCoverage, SEXP method, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP gaps(SEXP x, SEXP t);

SEXP firstSeqsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP first_x, SEXP first_y);

SEXP firstSeqsGapsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP t, SEXP first_x, SEXP first_y);

SEXP firstSeqsPosEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP t, SEXP first_x, SEXP first_y);

SEXP similarities(SEXP res, SEXP widths1, SEXP widths2, SEXP terminalGaps, SEXP penalizeGapLetters, SEXP minCoverage, SEXP method, SEXP nThreads);

SEXP overlap(SEXP res, SEXP widths1, SEXP widths2);

SEXP cophenetic(SEXP Index1, SEXP N, SEXP D, SEXP H);

// Cluster.c

SEXP cluster(SEXP x, SEXP cutoff, SEXP method, SEXP l, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP reclusterUPGMA(SEXP x, SEXP cutoff);

SEXP clusterNJ(SEXP x, SEXP cutoff, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP reclusterNJ(SEXP x, SEXP cutoff);

SEXP adjustHeights(SEXP x);

// ClusterML.c

SEXP clusterML(SEXP x, SEXP y, SEXP model, SEXP branches, SEXP lengths, SEXP states, SEXP type, SEXP weights, SEXP uLengths, SEXP nThreads);

SEXP expM(SEXP x, SEXP model, SEXP type);

// DesignProbes.c

SEXP designProbes(SEXP x, SEXP max_pl, SEXP min_pl, SEXP max_c, SEXP numMMs, SEXP numPs, SEXP st, SEXP en, SEXP max_ov, SEXP h_percent, SEXP min_f, SEXP max_f, SEXP minS, SEXP verbose, SEXP pBar, SEXP nThreads);

// CommonGaps.c

SEXP commonGaps(SEXP x);

// Utils.c

SEXP multiMatch(SEXP x, SEXP y, SEXP z);

SEXP multiMatchUpper(SEXP x, SEXP y, SEXP z);

SEXP multiMatchCharNotNA(SEXP x);

SEXP intMatch(SEXP x, SEXP y);

SEXP matchOrder(SEXP x, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP matchRanges(SEXP x, SEXP y, SEXP wordSize, SEXP maxLength, SEXP threshold);

SEXP boundedMatches(SEXP x, SEXP bl, SEXP bu);

SEXP intMatchOnce(SEXP x, SEXP y, SEXP o1, SEXP o2);

SEXP intMatchSelfOnce(SEXP x, SEXP o1);

SEXP matchListsDual(SEXP x, SEXP y, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP groupMax(SEXP x, SEXP y, SEXP z);

SEXP matchOverlap(SEXP x, SEXP y, SEXP v, SEXP wordSize, SEXP nThreads);

SEXP countHits(SEXP x, SEXP v, SEXP nThreads);

SEXP sumBins(SEXP v, SEXP bins);

SEXP xorShift(SEXP seed, SEXP base);

SEXP selectIndices(SEXP P, SEXP select, SEXP initial, SEXP final, SEXP ordering, SEXP num);

SEXP selectGroups(SEXP ordering, SEXP initial, SEXP final, SEXP num);

SEXP firstHit(SEXP x, SEXP y);

SEXP maxGroup(SEXP groups, SEXP weights);

SEXP sortedUnique(SEXP v);

SEXP splitPartitions(SEXP order, SEXP partition, SEXP var, SEXP minSize, SEXP minSplit);

// ReplaceChars.c

SEXP replaceChars(SEXP x, SEXP r, SEXP t);

SEXP replaceChar(SEXP x, SEXP c, SEXP r);

SEXP replaceGaps(SEXP x, SEXP y, SEXP start, SEXP type);

// TerminalMismatch.c

SEXP terminalMismatch(SEXP p, SEXP t, SEXP cutoff, SEXP mGaps, SEXP nThreads);

// NNLS.c

SEXP NNLS(SEXP row, SEXP col, SEXP value, SEXP nrows, SEXP ncols, SEXP b, SEXP tol, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP sparseMult(SEXP row, SEXP col, SEXP value, SEXP nrows, SEXP ncols, SEXP b);

// CalculateDeltaG.c

SEXP calculateDeltaG(SEXP p, SEXP t, SEXP deltaGrules);

SEXP calculateHairpinDeltaG(SEXP x, SEXP arms, SEXP deltaGrules);

// CalculateFISH.c

SEXP calculateFISH(SEXP probes, SEXP targets);

// AlignProfiles.c

SEXP alignProfiles(SEXP p, SEXP s, SEXP type, SEXP subMatrix, SEXP dbnMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP norm, SEXP nThreads);

SEXP alignProfilesAA(SEXP p, SEXP s, SEXP subMatrix, SEXP hecMatrix, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP norm, SEXP nThreads);

// EnumerateSequence.c

SEXP enumerateSequence(SEXP x, SEXP wordSize, SEXP mask, SEXP nThreads);

SEXP enumerateSequenceAA(SEXP x, SEXP wordSize);

SEXP enumerateGappedSequence(SEXP x, SEXP wordSize, SEXP ordering);

SEXP enumerateGappedSequenceAA(SEXP x, SEXP wordSize, SEXP ordering);

SEXP enumerateSequenceReducedAA(SEXP x, SEXP wordSize, SEXP alphabet, SEXP mask, SEXP nThreads);

SEXP alphabetSizeReducedAA(SEXP x, SEXP alphabet);

SEXP alphabetSize(SEXP x);

// Compositions.c

SEXP composition(SEXP x);

SEXP positionWeightMatrix(SEXP x, SEXP begins, SEXP ends, SEXP width);

// IntDist.c

SEXP intDist(SEXP x, SEXP levels, SEXP bins, SEXP maxBins, SEXP numRows, SEXP totRows, SEXP power);

// MeltPolymer.c

SEXP meltPolymer(SEXP x, SEXP temps, SEXP ions, SEXP output);

// InsertGaps.c

SEXP insertGaps(SEXP x, SEXP positions, SEXP lengths, SEXP type, SEXP nThreads);

// ExpandAmbiguities.c

SEXP expandAmbiguities(SEXP x, SEXP c);

// RemoveGaps.c

SEXP removeCommonGaps(SEXP x, SEXP type, SEXP mask, SEXP nThreads);

SEXP removeGaps(SEXP x, SEXP type, SEXP mask, SEXP nThreads);

// PredictHEC.c

SEXP predictHEC(SEXP x, SEXP windowSize, SEXP background, SEXP HEC_MI1, SEXP HEC_MI2, SEXP output);

// AssignIndels.c

SEXP clearIns(SEXP x);

SEXP all(SEXP x);

SEXP any(SEXP x);

// ConsolidateGaps.c

SEXP consolidateGaps(SEXP x, SEXP type);

// FindFrameshifts.c

SEXP findFrameshifts(SEXP t, SEXP l, SEXP f, SEXP index, SEXP oindex, SEXP maxComp, SEXP go, SEXP ge, SEXP fs, SEXP minD, SEXP maxD, SEXP subMatrix, SEXP verbose, SEXP pBar);

// Order.c

SEXP radixOrder(SEXP x, SEXP ascending, SEXP keepNAs);

SEXP dereplicate(SEXP x, SEXP o);

// ChainSegments.c

SEXP fillOverlaps(SEXP m, SEXP n);

SEXP indexByContig(SEXP starts, SEXP ends, SEXP order, SEXP index, SEXP widths);

SEXP chainSegments(SEXP x_s, SEXP x_e, SEXP x_i, SEXP x_f, SEXP y_s, SEXP y_e, SEXP y_i, SEXP y_f, SEXP weights, SEXP sepCost, SEXP gapCost, SEXP shiftCost, SEXP codingCost, SEXP maxSep, SEXP maxGap, SEXP ordering, SEXP minScore, SEXP minW, SEXP allowOverlap);

SEXP extendSegments(SEXP X, SEXP W1, SEXP W2, SEXP S1, SEXP S2, SEXP O1P, SEXP O1N, SEXP O2P, SEXP O2N, SEXP S, SEXP maxDrop, SEXP INDEX1, SEXP INDEX2);

// Translate.c

SEXP basicTranslate(SEXP x, SEXP code, SEXP starts);

// Import.c

SEXP collapse(SEXP x, SEXP index1, SEXP index2);

SEXP extractFields(SEXP x, SEXP fields, SEXP starts, SEXP ends);

// Compression.c

SEXP nbit(SEXP x, SEXP y, SEXP compRepeats, SEXP nThreads);

SEXP qbit(SEXP x, SEXP y, SEXP nThreads);

SEXP decompress(SEXP x, SEXP nThreads);

// Diff.c

SEXP intDiff(SEXP x);

// MovingAverage.c

SEXP movAvg(SEXP x, SEXP type, SEXP alpha, SEXP thresh, SEXP maxAvg, SEXP start, SEXP end);

// GetPools.c

SEXP getPools(SEXP x);

// PredictDBN.c

SEXP predictDBN(SEXP x, SEXP output, SEXP minOccupancy, SEXP impact, SEXP avgProdCorr, SEXP slope, SEXP shift, SEXP weights, SEXP pseudoknots, SEXP threshold, SEXP patterns, SEXP verbose, SEXP pBar, SEXP nThreads);

// InformationContent.c

SEXP informationContent(SEXP p, SEXP nS, SEXP correction, SEXP randomBackground);

SEXP informationContentAA(SEXP p, SEXP nS, SEXP correction, SEXP randomBackground);

// VectorSums.c

SEXP vectorSum(SEXP x, SEXP y, SEXP z, SEXP b);

SEXP parallelMatch(SEXP x, SEXP y, SEXP indices, SEXP a, SEXP b, SEXP pos, SEXP rng, SEXP nThreads);

// GeneFinding.c

SEXP getORFs(SEXP x, SEXP start_codons, SEXP stop_codons, SEXP min_gene_length, SEXP allow_edges);

SEXP codonModel(SEXP x, SEXP orftable, SEXP stop_codons, SEXP min_orf_length, SEXP coding_scores);

SEXP scoreCodonModel(SEXP x, SEXP orftable, SEXP codon_scores);

SEXP dicodonModel(SEXP x, SEXP orftable, SEXP stop_codons);

SEXP startCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP start_codons);

SEXP scoreStartCodonModel(SEXP x, SEXP orftable, SEXP start_scores);

SEXP initialCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP initial_codons);

SEXP scoreInitialCodonModel(SEXP x, SEXP orftable, SEXP ini_scores);

SEXP getRegion(SEXP x, SEXP orftable, SEXP width, SEXP offset, SEXP toStart);

SEXP autocorrelationModel(SEXP x, SEXP orftable, SEXP indices, SEXP aatable);

SEXP scoreAutocorrelationModel(SEXP x, SEXP orftable, SEXP codon_scores, SEXP aatable);

SEXP nucleotideBiasModel(SEXP x, SEXP orftable, SEXP indices, SEXP positions);

SEXP scoreNucleotideBiasModel(SEXP x, SEXP orftable, SEXP nuc_scores);

SEXP upstreamMotifModel(SEXP x, SEXP orftable, SEXP indices, SEXP begin, SEXP distance, SEXP size);

SEXP scoreUpstreamMotifModel(SEXP x, SEXP orftable, SEXP motif_scores, SEXP begin, SEXP distance, SEXP size);

SEXP runLengthModel(SEXP x, SEXP orftable, SEXP codon_scores);

SEXP scoreRunLengthModel(SEXP x, SEXP orftable, SEXP codon_scores, SEXP run_scores);

SEXP stopCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP stop_codons);

SEXP scoreStopCodonModel(SEXP x, SEXP orftable, SEXP stop_scores);

SEXP terminationCodonModel(SEXP x, SEXP orftable, SEXP indices, SEXP terminal_codons);

SEXP scoreTerminationCodonModel(SEXP x, SEXP orftable, SEXP ter_scores);

SEXP codonFrequencies(SEXP x, SEXP orftable, SEXP indices);

SEXP unicodonModel(SEXP x, SEXP orftable, SEXP stop_codons);

SEXP chainGenes(SEXP orftable, SEXP topScore, SEXP topLength, SEXP scoreIntergenic, SEXP maxOverlapSame, SEXP maxOverlapOpposite, SEXP maxFracOverlap, SEXP sameScores, SEXP oppoScores);

SEXP longestORFs(SEXP orftable);

SEXP getIndex(SEXP start1, SEXP start2, SEXP len, SEXP score);

SEXP getBounds(SEXP widths, SEXP start, SEXP end, SEXP minL, SEXP maxL, SEXP lenScores, SEXP kmer, SEXP Ksize, SEXP negOk, SEXP minS, SEXP partS, SEXP minC);

SEXP addIfElse(SEXP vec, SEXP index, SEXP scores);

SEXP kmerScores(SEXP oligos, SEXP ints, SEXP windowSize, SEXP kSize);

SEXP getHits(SEXP starts, SEXP ends, SEXP left1, SEXP left2, SEXP right1, SEXP right2, SEXP deltaG);

SEXP couplingModel(SEXP x, SEXP orftable, SEXP indices, SEXP aatable, SEXP maxDist);

SEXP scoreCouplingModel(SEXP x, SEXP orftable, SEXP coupling_scores, SEXP aatable);

SEXP maxPerORF(SEXP orftable, SEXP scores);

SEXP scorePWM(SEXP pwm, SEXP x, SEXP minScore, SEXP nThreads);

SEXP scoreTopPWM(SEXP pwm, SEXP x, SEXP begin, SEXP positions, SEXP nThreads);

SEXP dist(SEXP x, SEXP nThreads);

// ClusterMP.c

SEXP clusterMP(SEXP x, SEXP z, SEXP s, SEXP sizes, SEXP scoreOnly, SEXP add, SEXP weights, SEXP nThreads);

// PairwiseAlignment.c

SEXP alignPair(SEXP x, SEXP y, SEXP s1, SEXP e1, SEXP s2, SEXP e2, SEXP go, SEXP ge, SEXP tg, SEXP maxLength, SEXP type, SEXP subMatrix, SEXP nThreads);
