# annogene

## Example
```
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

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func usage() {
	fmt.Printf("\nProgram: BiTransformFastQ (Tools for FASTQ C to T or G to A Transform)\nVersion: 0.1.1-20170630\n\nUsage:\tBiTransformFastQ -inFQ sample1_P_R1.fq.gz -TF C2T -o /abspath/outdir\n\n\tThe out file is /abspath/outdir/sample1_P_R1.fq.gz_C2T.fq.gz\n")
	fmt.Printf("Command:\n")

	fmt.Printf("    -inFQ          faseq.gz\n")
	fmt.Printf("    -TF            C2T or G2A\n")
	fmt.Printf("    -o             outdir\n")
	os.Exit(1)
}

func main() {
	infq := flag.String("inFQ", "", "fastq.gz")
	transf := flag.String("TF", "", "C2T / G2A")
	outdir := flag.String("o", "", "outdir")
	flag.Parse()
	if *infq == "" || *transf == "" || *outdir == "" {
		usage()
	}

	file, err := os.Open(*infq)
	check(err)
	gz, err := gzip.NewReader(file)
	check(err)

	defer file.Close()
	defer gz.Close()

	r := fastq.NewReader(gz)
	sc := fastq.NewScanner(r)

	outfqgz := fmt.Sprintf("%s/%s_%s.fq.gz", *outdir, filepath.Base(*infq), *transf)

	fo, err := os.OpenFile(outfqgz, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0660)
	check(err)
	ogz := gzip.NewWriter(fo)
	check(err)

	defer fo.Close()
	defer ogz.Close()

	w := fastq.NewWriter(ogz)

	switch *transf {
	case "C2T":
		for sc.Next() {
			CTreads := fastq.C2T(sc.Seq())
			_, eer := w.Write(CTreads)
			check(eer)
		}
		if err := sc.Error(); err != nil {
			log.Fatalf("failed to read fastq: %v", err)
		}
	case "G2A":
		for sc.Next() {
			GAreads := fastq.G2A(sc.Seq())
			_, eer := w.Write(GAreads)
			check(eer)
		}
		if err := sc.Error(); err != nil {
			log.Fatalf("failed to read fastq: %v", err)
		}
	default:
		usage()
	}
}
```


```
package main

import (
	"flag"
	"fmt"
	"github.com/seqyuan/annogene/io/fastq"
	"log"
	"os"
	)

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func usage() {
	fmt.Printf("\nProgram: cut fastq length\n")
	fmt.Printf("Command:\n")
	fmt.Printf("    -inFQ          in.faseq\n")
	fmt.Printf("    -c             cut length from 5'\n")
	fmt.Printf("    -o             outfile.fastq\n")
	os.Exit(1)
}

func main() {
	infq := flag.String("inFQ", "", "test.fastq")
	cutLen := flag.Int("c", 30, "30 or other int umber")
	outfile := flag.String("o", "", "outfile.fastq")
	flag.Parse()
	if *infq == "" || *outfile == "" {
		usage()
	}

	file, err := os.Open(*infq)
	check(err)

	defer file.Close()

	r := fastq.NewReader(file)
	sc := fastq.NewScanner(r)

	fo, err := os.OpenFile(*outfile, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0660)
	check(err)

	defer fo.Close()

	w := fastq.NewWriter(fo)

	for sc.Next() {
		CuTreads := fastq.CutLen(sc.Seq(), *cutLen)
		_, eer := w.Write(CuTreads)
		check(eer)
		}
	if err := sc.Error(); err != nil {
		log.Fatalf("failed to read fastq: %v", err)
	}

}

```
