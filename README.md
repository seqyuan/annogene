# annogene

##Example
```
package main

import (
	"compress/gzip"
	//	"fmt"
	"github.com/seqyuan/annogene/io/fastq"
	"log"
	"os"
)

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func main() {
	filename := "D:\\GoWorkspace\\src\\KG1a-o-twist_P_R1.fq.gz"
	file, err := os.Open(filename)
	check(err)

	gz, err := gzip.NewReader(file)
	check(err)

	defer file.Close()
	defer gz.Close()

	r := fastq.NewReader(gz)
	sc := fastq.NewScanner(r)

	outfile := "D:\\GoWorkspace\\src\\KG1a-o-twist_P_R1.fq.gz_C_to_T.fastq_1.gz"
	fo, err := os.OpenFile(outfile, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0660)
	check(err)
	ogz := gzip.NewWriter(fo)
	check(err)

	defer fo.Close()
	defer ogz.Close()

	w := fastq.NewWriter(ogz)

	for sc.Next() {
		GA := fastq.G2A(sc.Seq())
		CTreads := fastq.C2T(GA)

		_, eer := w.Write(CTreads)
		check(eer)
	}

	if err := sc.Error(); err != nil {
		log.Fatalf("failed to read fastq: %v", err)
	}

}
```