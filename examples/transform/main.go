// Package main demonstrates how to use the annogene package for FASTQ transformation
package main

import (
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/seqyuan/annogene/io/fastq"
)

func main() {
	// Parse command line flags
	infq := flag.String("inFQ", "", "Input FASTQ file (supports .gz)")
	transf := flag.String("TF", "", "Transformation: C2T or G2A")
	outdir := flag.String("o", "", "Output directory")
	flag.Parse()

	// Validate required parameters
	if *infq == "" || *transf == "" || *outdir == "" {
		fmt.Println("Usage: transform -inFQ <input.fastq> -TF <C2T|G2A> -o <output_dir>")
		flag.PrintDefaults()
		os.Exit(1)
	}

	// Validate transformation type
	if *transf != "C2T" && *transf != "G2A" {
		log.Fatalf("Invalid transformation: %s. Must be C2T or G2A", *transf)
	}

	// Open input file
	file, err := os.Open(*infq)
	if err != nil {
		log.Fatalf("Failed to open input file: %v", err)
	}
	defer file.Close()

	// Handle gzipped files
	var reader io.Reader = file
	if filepath.Ext(*infq) == ".gz" {
		gz, err := gzip.NewReader(file)
		if err != nil {
			log.Fatalf("Failed to create gzip reader: %v", err)
		}
		defer gz.Close()
		reader = gz
	}

	// Create FASTQ reader and scanner
	r := fastq.NewReader(reader)
	scanner := fastq.NewScanner(r)

	// Prepare output file path
	baseName := filepath.Base(*infq)
	if filepath.Ext(baseName) == ".gz" {
		baseName = strings.TrimSuffix(baseName, ".gz")
	}
	outPath := filepath.Join(*outdir, fmt.Sprintf("%s_%s.fastq", baseName, *transf))

	// Create output directory if it doesn't exist
	if err := os.MkdirAll(*outdir, 0755); err != nil {
		log.Fatalf("Failed to create output directory: %v", err)
	}

	// Create output file
	outFile, err := os.Create(outPath)
	if err != nil {
		log.Fatalf("Failed to create output file: %v", err)
	}
	defer outFile.Close()

	writer := fastq.NewWriter(outFile)

	// Process sequences
	var processedCount int
	switch *transf {
	case "C2T":
		log.Printf("Converting C to T in sequences...")
		for scanner.Next() {
			converted := fastq.C2T(scanner.Seq())
			if _, err := writer.Write(converted); err != nil {
				log.Printf("Warning: Failed to write sequence: %v", err)
				continue
			}
			processedCount++
		}
	case "G2A":
		log.Printf("Converting G to A in sequences...")
		for scanner.Next() {
			converted := fastq.G2A(scanner.Seq())
			if _, err := writer.Write(converted); err != nil {
				log.Printf("Warning: Failed to write sequence: %v", err)
				continue
			}
			processedCount++
		}
	}

	// Check for scanner errors
	if err := scanner.Error(); err != nil {
		log.Fatalf("Failed to read FASTQ: %v", err)
	}

	log.Printf("Successfully processed %d sequences", processedCount)
	log.Printf("Output written to: %s", outPath)
}
