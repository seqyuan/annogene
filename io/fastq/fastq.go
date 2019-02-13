package fastq

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"log"
	"fmt"
)

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

type sReader struct {
	r *bufio.Reader
}

func NewReader(r io.Reader) *sReader {
	return &sReader{
		r: bufio.NewReader(r),
	}
}

type Reader interface {
	// Read reads a seq.Sequence, returning the sequence and any error that
	// occurred during the read.
	Read() (Sequence, error)
}

type Scanner struct {
	r   Reader
	seq Sequence
	err error
}

func (s *Scanner) Next() bool {
	if s.err != nil {
		return false
	}
	s.seq, s.err = s.r.Read()
	return s.err == nil
}

// Error returns the first non-EOF error that was encountered by the Scanner.
func (s *Scanner) Error() error {
	if s.err == io.EOF {
		return nil
	}
	return s.err
}

// Seq returns the most recent sequence read by a call to Next.
func (s *Scanner) Seq() Sequence { return s.seq }

// NewScanner returns a Scanner to read from r.
func NewScanner(r Reader) *Scanner { return &Scanner{r: r} }

type Sequence struct {
	Id1     []byte
	Letters []byte
	Id2     []byte
	Quality []byte
}

func (reads *Sequence) SetId1(id []byte) error          { reads.Id1 = id; return nil }
func (reads *Sequence) SetLetters(letters []byte) error { reads.Letters = letters; return nil }
func (reads *Sequence) SetId2(id2 []byte) error         { reads.Id2 = id2; return nil }
func (reads *Sequence) SetQuality(quality []byte) error { reads.Quality = quality; return nil }

func C2T(reads Sequence) Sequence {
	cc := []byte("C")
	tt := []byte("T")
	ss2 := bytes.Replace(reads.Letters, cc, tt, -1)
	reads.Letters = ss2
	return reads
}

func G2A(reads Sequence) Sequence {
	gg := []byte("G")
	aa := []byte("A")
	ss2 := bytes.Replace(reads.Letters, gg, aa, -1)
	reads.Letters = ss2
	return reads
}

func CutLen(reads Sequence, leng int)( reads2 Sequence){
	reads2.Id1 = reads.Id1
	if len(reads.Letters) > leng{
		reads2.Letters = reads.Letters[0:leng]
		reads2.Quality = reads.Quality[0:leng]
	}else{
		reads2.Letters = reads.Letters
		reads2.Quality = reads.Quality
	}
	reads2.Id2 = reads.Id2
	return 
}


func (r *sReader) Read() (Sequence, error) {
	const (
		id1 = iota
		letters
		id2
		quality
	)

	var (
		buff, line []byte
		isPrefix   bool
		state      int
		err        error
		reads      Sequence
	)
	fmt.Println(001,string(reads.Id1))
loop:
	for {
		fmt.Println(111,string(reads.Id1))
		buff, isPrefix, err = r.r.ReadLine()
		if err != nil {
			if state == quality && err == io.EOF {
				err = nil
				break
			}
			return reads, err
		}
		line = append(line, buff...)
		if isPrefix {
			continue
		}

		line = bytes.TrimSpace(line)
		switch {
		case state == id1 && maybeID1(line):
			state = letters
			err = reads.SetId1(line)
			check(err)

		case state == id2 && maybeID2(line):
			state = quality
			ii := append([]byte(nil), line...)
			err = reads.SetId2(ii)
			check(err)
			fmt.Println(113,string(reads.Id1))

		case state == letters && len(line) > 0:
			fmt.Println(112,string(reads.Id1))
			if maybeID2(line) && (len(line) == 1 || bytes.Compare(reads.Id1[1:], line[1:]) == 0) {
				state = quality
				break
			}
			ll := append([]byte(nil), line...)
			err = reads.SetLetters(ll)
			check(err)

			state = id2
		case state == quality:
			if len(line) == 0 && len(reads.Letters) != 0 {
				continue
			}
			break loop
		}
		line = line[:0]
	}

	line = bytes.Join(bytes.Fields(line), nil)
	if len(line) != len(reads.Letters) {
		return reads, errors.New("fastq: sequence/quality length mismatch")
	}
	err = reads.SetQuality(line)
	check(err)
	fmt.Println(114,string(reads.Id1))
	return reads, err
}

func maybeID1(l []byte) bool { return len(l) > 0 && l[0] == '@' }
func maybeID2(l []byte) bool { return len(l) > 0 && l[0] == '+' }

// Fastq sequence format writer type.
type Writer struct {
	w io.Writer
}

// Returns a new fastq format writer using w.
func NewWriter(w io.Writer) *Writer {
	return &Writer{
		w: w,
	}
}

// Write a single sequence and return the number of bytes written and any error.
func (w *Writer) Write(s Sequence) (n int, err error) {
	var _n int

	n, err = w.w.Write(s.Id1)
	if n += _n; err != nil {
		return
	}
	_n, err = w.w.Write([]byte{'\n'})
	if n += _n; err != nil {
		return
	}

	n, err = w.w.Write(s.Letters)
	if n += _n; err != nil {
		return
	}
	_n, err = w.w.Write([]byte{'\n'})
	if n += _n; err != nil {
		return
	}

	n, err = w.w.Write(s.Id2)
	if n += _n; err != nil {
		return
	}
	_n, err = w.w.Write([]byte{'\n'})
	if n += _n; err != nil {
		return
	}

	n, err = w.w.Write(s.Quality)
	if n += _n; err != nil {
		return
	}
	_n, err = w.w.Write([]byte{'\n'})
	if n += _n; err != nil {
		return
	}

	return
}
