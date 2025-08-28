// Package fastq provides functionality for reading and writing FASTQ format files.
// FASTQ is a text-based format for storing biological sequences and their quality scores.
package fastq

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Common errors
var (
	ErrSequenceQualityMismatch = errors.New("fastq: sequence/quality length mismatch")
	ErrInvalidRegionFormat     = errors.New("fastq: invalid region format")
)

// Reader interface defines methods for reading FASTQ sequences
type Reader interface {
	// Read reads a Sequence, returning the sequence and any error that
	// occurred during the read.
	Read() (Sequence, error)
}

// Writer interface defines methods for writing FASTQ sequences
type Writer interface {
	// Write writes a single sequence and returns the number of bytes written and any error.
	Write(s Sequence) (n int, err error)
}

// Sequence represents a single FASTQ sequence with ID, sequence letters, and quality scores
type Sequence struct {
	ID1     []byte // First identifier line (starts with @)
	Letters []byte // Sequence letters
	ID2     []byte // Second identifier line (starts with +)
	Quality []byte // Quality scores
}

// SetID1 sets the first identifier for the sequence
func (s *Sequence) SetID1(id []byte) error {
	if id == nil {
		return errors.New("id cannot be nil")
	}
	s.ID1 = make([]byte, len(id))
	copy(s.ID1, id)
	return nil
}

// SetLetters sets the sequence letters
func (s *Sequence) SetLetters(letters []byte) error {
	if letters == nil {
		return errors.New("letters cannot be nil")
	}
	s.Letters = make([]byte, len(letters))
	copy(s.Letters, letters)
	return nil
}

// SetID2 sets the second identifier for the sequence
func (s *Sequence) SetID2(id2 []byte) error {
	if id2 == nil {
		return errors.New("id2 cannot be nil")
	}
	s.ID2 = make([]byte, len(id2))
	copy(s.ID2, id2)
	return nil
}

// SetQuality sets the quality scores for the sequence
func (s *Sequence) SetQuality(quality []byte) error {
	if quality == nil {
		return errors.New("quality cannot be nil")
	}
	s.Quality = make([]byte, len(quality))
	copy(s.Quality, quality)
	return nil
}

// C2T converts all C bases to T bases in the sequence
func C2T(seq Sequence) Sequence {
	result := Sequence{}
	result.ID1 = make([]byte, len(seq.ID1))
	copy(result.ID1, seq.ID1)
	result.ID2 = make([]byte, len(seq.ID2))
	copy(result.ID2, seq.ID2)
	
	// Convert C to T
	result.Letters = bytes.ReplaceAll(seq.Letters, []byte("C"), []byte("T"))
	result.Quality = make([]byte, len(seq.Quality))
	copy(result.Quality, seq.Quality)
	
	return result
}

// G2A converts all G bases to A bases in the sequence
func G2A(seq Sequence) Sequence {
	result := Sequence{}
	result.ID1 = make([]byte, len(seq.ID1))
	copy(result.ID1, seq.ID1)
	result.ID2 = make([]byte, len(seq.ID2))
	copy(result.ID2, seq.ID2)
	
	// Convert G to A
	result.Letters = bytes.ReplaceAll(seq.Letters, []byte("G"), []byte("A"))
	result.Quality = make([]byte, len(seq.Quality))
	copy(result.Quality, seq.Quality)
	
	return result
}

// CutLen truncates the sequence to the specified length from the 5' end
func CutLen(seq Sequence, length int) Sequence {
	if length <= 0 {
		return seq
	}
	
	result := Sequence{}
	result.ID1 = make([]byte, len(seq.ID1))
	copy(result.ID1, seq.ID1)
	result.ID2 = make([]byte, len(seq.ID2))
	copy(result.ID2, seq.ID2)
	
	if len(seq.Letters) > length {
		result.Letters = seq.Letters[:length]
		result.Quality = seq.Quality[:length]
	} else {
		result.Letters = make([]byte, len(seq.Letters))
		copy(result.Letters, seq.Letters)
		result.Quality = make([]byte, len(seq.Quality))
		copy(result.Quality, seq.Quality)
	}
	
	return result
}

// ExtractRegion extracts specified regions from the sequence
func ExtractRegion(seq Sequence, regions string) (Sequence, error) {
	if regions == "" {
		return seq, nil
	}
	
	result := Sequence{}
	result.ID1 = make([]byte, len(seq.ID1))
	copy(result.ID1, seq.ID1)
	result.ID2 = make([]byte, len(seq.ID2))
	copy(result.ID2, seq.ID2)
	
	var extractedLetters []byte
	var extractedQuality []byte
	
	// Split regions string to get each region's start and end indices
	regionPairs := strings.Split(regions, ",")
	for _, pair := range regionPairs {
		// Split each region's start and end indices
		indexRange := strings.Split(pair, ":")
		if len(indexRange) != 2 {
			return result, fmt.Errorf("%w: %s", ErrInvalidRegionFormat, pair)
		}
		
		start, err := strconv.Atoi(strings.TrimSpace(indexRange[0]))
		if err != nil {
			return result, fmt.Errorf("invalid start index: %w", err)
		}
		
		end, err := strconv.Atoi(strings.TrimSpace(indexRange[1]))
		if err != nil {
			return result, fmt.Errorf("invalid end index: %w", err)
		}
		
		// Validate indices
		if start < 0 || end > len(seq.Letters) || start >= end {
			return result, fmt.Errorf("invalid region indices: start=%d, end=%d, sequence_length=%d", start, end, len(seq.Letters))
		}
		
		// Extract substrings and add to slices
		extractedLetters = append(extractedLetters, seq.Letters[start:end]...)
		extractedQuality = append(extractedQuality, seq.Quality[start:end]...)
	}
	
	result.Letters = extractedLetters
	result.Quality = extractedQuality
	
	return result, nil
}

// Scanner provides a convenient interface for reading FASTQ sequences
type Scanner struct {
	r   Reader
	seq Sequence
	err error
}

// NewScanner returns a new Scanner to read from r
func NewScanner(r Reader) *Scanner {
	return &Scanner{r: r}
}

// Next advances the Scanner to the next sequence
func (s *Scanner) Next() bool {
	if s.err != nil {
		return false
	}
	s.seq, s.err = s.r.Read()
	return s.err == nil
}

// Error returns the first non-EOF error that was encountered by the Scanner
func (s *Scanner) Error() error {
	if s.err == io.EOF {
		return nil
	}
	return s.err
}

// Seq returns the most recent sequence read by a call to Next
func (s *Scanner) Seq() Sequence { return s.seq }

// fastqReader implements the Reader interface for FASTQ files
type fastqReader struct {
	r *bufio.Reader
}

// NewReader returns a new FASTQ reader
func NewReader(r io.Reader) Reader {
	return &fastqReader{
		r: bufio.NewReader(r),
	}
}

// Read reads a single FASTQ sequence
func (r *fastqReader) Read() (Sequence, error) {
	const (
		stateID1    = iota
		stateLetters
		stateID2
		stateQuality
	)

	var (
		buff, line []byte
		isPrefix   bool
		state      int
		err        error
		seq        Sequence
	)

loop:
	for {
		buff, isPrefix, err = r.r.ReadLine()
		if err != nil {
			if state == stateQuality && err == io.EOF {
				err = nil
				break
			}
			return seq, err
		}
		line = append(line, buff...)
		if isPrefix {
			continue
		}

		line = bytes.TrimSpace(line)
		switch {
		case state == stateID1 && maybeID1(line):
			state = stateLetters
			if err := seq.SetID1(line); err != nil {
				return seq, fmt.Errorf("failed to set ID1: %w", err)
			}

		case state == stateID2 && maybeID2(line):
			state = stateQuality
			if err := seq.SetID2(line); err != nil {
				return seq, fmt.Errorf("failed to set ID2: %w", err)
			}

		case state == stateLetters && len(line) > 0:
			if maybeID2(line) && (len(line) == 1 || bytes.Equal(seq.ID1[1:], line[1:])) {
				state = stateQuality
				break
			}
			if err := seq.SetLetters(line); err != nil {
				return seq, fmt.Errorf("failed to set letters: %w", err)
			}
			state = stateID2
			
		case state == stateQuality:
			if len(line) == 0 && len(seq.Letters) != 0 {
				continue
			}
			break loop
		}
		line = line[:0]
	}

	line = bytes.Join(bytes.Fields(line), nil)
	if len(line) != len(seq.Letters) {
		return seq, ErrSequenceQualityMismatch
	}
	
	if err := seq.SetQuality(line); err != nil {
		return seq, fmt.Errorf("failed to set quality: %w", err)
	}
	
	return seq, nil
}

// maybeID1 checks if a line might be an ID1 line (starts with @)
func maybeID1(l []byte) bool { return len(l) > 0 && l[0] == '@' }

// maybeID2 checks if a line might be an ID2 line (starts with +)
func maybeID2(l []byte) bool { return len(l) > 0 && l[0] == '+' }

// fastqWriter implements the Writer interface for FASTQ files
type fastqWriter struct {
	w io.Writer
}

// NewWriter returns a new FASTQ writer
func NewWriter(w io.Writer) Writer {
	return &fastqWriter{w: w}
}

// Write writes a single sequence and returns the number of bytes written and any error
func (w *fastqWriter) Write(s Sequence) (n int, err error) {
	var written int

	// Write ID1
	written, err = w.w.Write(s.ID1)
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write ID1: %w", err)
	}
	
	// Write newline
	written, err = w.w.Write([]byte{'\n'})
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write newline: %w", err)
	}

	// Write sequence letters
	written, err = w.w.Write(s.Letters)
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write letters: %w", err)
	}
	
	// Write newline
	written, err = w.w.Write([]byte{'\n'})
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write newline: %w", err)
	}

	// Write ID2
	written, err = w.w.Write(s.ID2)
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write ID2: %w", err)
	}
	
	// Write newline
	written, err = w.w.Write([]byte{'\n'})
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write newline: %w", err)
	}

	// Write quality scores
	written, err = w.w.Write(s.Quality)
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write quality: %w", err)
	}
	
	// Write final newline
	written, err = w.w.Write([]byte{'\n'})
	n += written
	if err != nil {
		return n, fmt.Errorf("failed to write newline: %w", err)
	}

	return n, nil
}
