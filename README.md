# Annogene

[![Go Version](https://img.shields.io/github/go-mod/go-version/seqyuan/annogene)](https://golang.org)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Go Report Card](https://goreportcard.com/badge/github.com/seqyuan/annogene)](https://goreportcard.com/report/github.com/seqyuan/annogene)

Annogene is a Go package for processing and transforming FASTQ format files. It provides functionality for base conversion, sequence length cutting, and region extraction from biological sequences.

## Features

- **FASTQ Reading/Writing**: Efficient parsing and writing of FASTQ format files
- **Base Conversion**: Convert C to T or G to A bases in sequences
- **Sequence Truncation**: Cut sequences to specified lengths from the 5' end
- **Region Extraction**: Extract specific regions from sequences
- **Memory Efficient**: Uses byte slices for optimal performance
- **Error Handling**: Comprehensive error handling with detailed error messages

## Installation

```bash
go get github.com/seqyuan/annogene
```

## Quick Start

```go
package main

import (
    "os"
    "github.com/seqyuan/annogene/io/fastq"
)

func main() {
    // Open a FASTQ file
    file, err := os.Open("sample.fastq")
    if err != nil {
        panic(err)
    }
    defer file.Close()

    // Create a reader
    reader := fastq.NewReader(file)
    scanner := fastq.NewScanner(reader)

    // Process sequences
    for scanner.Next() {
        seq := scanner.Seq()
        
        // Convert C to T
        converted := fastq.C2T(seq)
        
        // Process the converted sequence...
    }

    if err := scanner.Error(); err != nil {
        panic(err)
    }
}
```

## API Reference

### Core Types

#### Sequence
```go
type Sequence struct {
    ID1     []byte // First identifier line (starts with @)
    Letters []byte // Sequence letters
    ID2     []byte // Second identifier line (starts with +)
    Quality []byte // Quality scores
}
```

#### Reader Interface
```go
type Reader interface {
    Read() (Sequence, error)
}
```

#### Writer Interface
```go
type Writer interface {
    Write(s Sequence) (n int, err error)
}
```

### Functions

#### Base Conversion
- `C2T(seq Sequence) Sequence` - Converts all C bases to T bases
- `G2A(seq Sequence) Sequence` - Converts all G bases to A bases

#### Sequence Manipulation
- `CutLen(seq Sequence, length int) Sequence` - Truncates sequence to specified length
- `ExtractRegion(seq Sequence, regions string) (Sequence, error)` - Extracts specified regions

#### Utility Functions
- `NewReader(r io.Reader) Reader` - Creates a new FASTQ reader
- `NewWriter(w io.Writer) Writer` - Creates a new FASTQ writer
- `NewScanner(r Reader) *Scanner` - Creates a new scanner for reading sequences

### Scanner
```go
type Scanner struct {
    // ... internal fields
}

func (s *Scanner) Next() bool      // Advance to next sequence
func (s *Scanner) Seq() Sequence   // Get current sequence
func (s *Scanner) Error() error    // Get any error encountered
```

## Examples

### Base Transformation Tool

```go
package main

import (
    "compress/gzip"
    "flag"
    "fmt"
    "github.com/seqyuan/annogene/io/fastq"
    "log"
    "os"
    "path/filepath"
)

func main() {
    infq := flag.String("inFQ", "", "Input FASTQ file")
    transf := flag.String("TF", "", "Transformation: C2T or G2A")
    outdir := flag.String("o", "", "Output directory")
    flag.Parse()

    if *infq == "" || *transf == "" || *outdir == "" {
        flag.Usage()
        os.Exit(1)
    }

    // Open input file
    file, err := os.Open(*infq)
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()

    // Handle gzipped files
    var reader io.Reader = file
    if filepath.Ext(*infq) == ".gz" {
        gz, err := gzip.NewReader(file)
        if err != nil {
            log.Fatal(err)
        }
        defer gz.Close()
        reader = gz
    }

    // Create FASTQ reader and scanner
    r := fastq.NewReader(reader)
    scanner := fastq.NewScanner(r)

    // Prepare output
    outPath := fmt.Sprintf("%s/%s_%s.fastq", *outdir, 
        filepath.Base(*infq), *transf)
    outFile, err := os.Create(outPath)
    if err != nil {
        log.Fatal(err)
    }
    defer outFile.Close()

    writer := fastq.NewWriter(outFile)

    // Process sequences
    switch *transf {
    case "C2T":
        for scanner.Next() {
            converted := fastq.C2T(scanner.Seq())
            if _, err := writer.Write(converted); err != nil {
                log.Printf("Failed to write sequence: %v", err)
            }
        }
    case "G2A":
        for scanner.Next() {
            converted := fastq.G2A(scanner.Seq())
            if _, err := writer.Write(converted); err != nil {
                log.Printf("Failed to write sequence: %v", err)
            }
        }
    default:
        log.Fatalf("Unknown transformation: %s", *transf)
    }

    if err := scanner.Error(); err != nil {
        log.Fatalf("Failed to read FASTQ: %v", err)
    }
}
```

### Sequence Length Cutting Tool

```go
package main

import (
    "flag"
    "github.com/seqyuan/annogene/io/fastq"
    "log"
    "os"
)

func main() {
    infq := flag.String("inFQ", "", "Input FASTQ file")
    cutLen := flag.Int("c", 30, "Cut length from 5' end")
    outfile := flag.String("o", "", "Output FASTQ file")
    flag.Parse()

    if *infq == "" || *outfile == "" {
        flag.Usage()
        os.Exit(1)
    }

    // Open input file
    file, err := os.Open(*infq)
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()

    // Create reader and scanner
    reader := fastq.NewReader(file)
    scanner := fastq.NewScanner(reader)

    // Create output file
    outFile, err := os.Create(*outfile)
    if err != nil {
        log.Fatal(err)
    }
    defer outFile.Close()

    writer := fastq.NewWriter(outFile)

    // Process sequences
    for scanner.Next() {
        cutSeq := fastq.CutLen(scanner.Seq(), *cutLen)
        if _, err := writer.Write(cutSeq); err != nil {
            log.Printf("Failed to write sequence: %v", err)
        }
    }

    if err := scanner.Error(); err != nil {
        log.Fatalf("Failed to read FASTQ: %v", err)
    }
}
```

## Development

### Prerequisites
- Go 1.21 or later

### Building
```bash
make build
```

### Testing
```bash
make test
```

### Code Formatting
```bash
make format
```

### Linting
```bash
make lint
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Version History

- **v0.0.2** - Current version with improved error handling and code structure
- **v0.0.1** - Initial release

## Acknowledgments

This package was developed for bioinformatics applications and sequence analysis workflows.
