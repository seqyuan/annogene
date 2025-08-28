package fastq

import (
	"bytes"
	"strings"
	"testing"
)

func TestSequence_SetID1(t *testing.T) {
	seq := &Sequence{}
	id := []byte("@test_id")
	
	err := seq.SetID1(id)
	if err != nil {
		t.Errorf("SetID1 failed: %v", err)
	}
	
	if !bytes.Equal(seq.ID1, id) {
		t.Errorf("Expected ID1 %v, got %v", id, seq.ID1)
	}
	
	// Test nil input
	err = seq.SetID1(nil)
	if err == nil {
		t.Error("Expected error for nil input")
	}
}

func TestSequence_SetLetters(t *testing.T) {
	seq := &Sequence{}
	letters := []byte("ACGT")
	
	err := seq.SetLetters(letters)
	if err != nil {
		t.Errorf("SetLetters failed: %v", err)
	}
	
	if !bytes.Equal(seq.Letters, letters) {
		t.Errorf("Expected Letters %v, got %v", letters, seq.Letters)
	}
	
	// Test nil input
	err = seq.SetLetters(nil)
	if err == nil {
		t.Error("Expected error for nil input")
	}
}

func TestSequence_SetID2(t *testing.T) {
	seq := &Sequence{}
	id2 := []byte("+test_id")
	
	err := seq.SetID2(id2)
	if err != nil {
		t.Errorf("SetID2 failed: %v", err)
	}
	
	if !bytes.Equal(seq.ID2, id2) {
		t.Errorf("Expected ID2 %v, got %v", id2, seq.ID2)
	}
	
	// Test nil input
	err = seq.SetID2(nil)
	if err == nil {
		t.Error("Expected error for nil input")
	}
}

func TestSequence_SetQuality(t *testing.T) {
	seq := &Sequence{}
	quality := []byte("!!!!")
	
	err := seq.SetQuality(quality)
	if err != nil {
		t.Errorf("SetQuality failed: %v", err)
	}
	
	if !bytes.Equal(seq.Quality, quality) {
		t.Errorf("Expected Quality %v, got %v", quality, seq.Quality)
	}
	
	// Test nil input
	err = seq.SetQuality(nil)
	if err == nil {
		t.Error("Expected error for nil input")
	}
}

func TestC2T(t *testing.T) {
	original := Sequence{
		ID1:     []byte("@test"),
		Letters: []byte("ACGT"),
		ID2:     []byte("+test"),
		Quality: []byte("!!!!"),
	}
	
	result := C2T(original)
	
	// Check that C was converted to T
	expected := []byte("ATGT")
	if !bytes.Equal(result.Letters, expected) {
		t.Errorf("Expected %v, got %v", expected, result.Letters)
	}
	
	// Check that other fields are preserved
	if !bytes.Equal(result.ID1, original.ID1) {
		t.Errorf("ID1 was modified")
	}
	if !bytes.Equal(result.ID2, original.ID2) {
		t.Errorf("ID2 was modified")
	}
	if !bytes.Equal(result.Quality, original.Quality) {
		t.Errorf("Quality was modified")
	}
	
	// Check that original is not modified
	if !bytes.Equal(original.Letters, []byte("ACGT")) {
		t.Errorf("Original sequence was modified")
	}
}

func TestG2A(t *testing.T) {
	original := Sequence{
		ID1:     []byte("@test"),
		Letters: []byte("ACGT"),
		ID2:     []byte("+test"),
		Quality: []byte("!!!!"),
	}
	
	result := G2A(original)
	
	// Check that G was converted to A
	expected := []byte("ACAT")
	if !bytes.Equal(result.Letters, expected) {
		t.Errorf("Expected %v, got %v", expected, result.Letters)
	}
	
	// Check that other fields are preserved
	if !bytes.Equal(result.ID1, original.ID1) {
		t.Errorf("ID1 was modified")
	}
	if !bytes.Equal(result.ID2, original.ID2) {
		t.Errorf("ID2 was modified")
	}
	if !bytes.Equal(result.Quality, original.Quality) {
		t.Errorf("Quality was modified")
	}
	
	// Check that original is not modified
	if !bytes.Equal(original.Letters, []byte("ACGT")) {
		t.Errorf("Original sequence was modified")
	}
}

func TestCutLen(t *testing.T) {
	original := Sequence{
		ID1:     []byte("@test"),
		Letters: []byte("ACGTACGT"),
		ID2:     []byte("+test"),
		Quality: []byte("!!!!!!!!"),
	}
	
	// Test cutting to shorter length
	result := CutLen(original, 4)
	expected := []byte("ACGT")
	if !bytes.Equal(result.Letters, expected) {
		t.Errorf("Expected %v, got %v", expected, result.Letters)
	}
	
	// Test cutting to longer length (should return original)
	result = CutLen(original, 10)
	if !bytes.Equal(result.Letters, original.Letters) {
		t.Errorf("Expected original length, got %v", result.Letters)
	}
	
	// Test cutting to zero or negative length
	result = CutLen(original, 0)
	if !bytes.Equal(result.Letters, original.Letters) {
		t.Errorf("Expected original for zero length")
	}
	
	result = CutLen(original, -1)
	if !bytes.Equal(result.Letters, original.Letters) {
		t.Errorf("Expected original for negative length")
	}
}

func TestExtractRegion(t *testing.T) {
	original := Sequence{
		ID1:     []byte("@test"),
		Letters: []byte("ACGTACGT"),
		ID2:     []byte("+test"),
		Quality: []byte("!!!!!!!!"),
	}
	
	// Test single region extraction
	result, err := ExtractRegion(original, "1:4")
	if err != nil {
		t.Errorf("ExtractRegion failed: %v", err)
	}
	
	expected := []byte("CGT")
	if !bytes.Equal(result.Letters, expected) {
		t.Errorf("Expected %v, got %v", expected, result.Letters)
	}
	
	// Test multiple region extraction
	result, err = ExtractRegion(original, "0:2,6:8")
	if err != nil {
		t.Errorf("ExtractRegion failed: %v", err)
	}
	
	expected = []byte("ACGT")
	if !bytes.Equal(result.Letters, expected) {
		t.Errorf("Expected %v, got %v", expected, result.Letters)
	}
	
	// Test empty regions
	result, err = ExtractRegion(original, "")
	if err != nil {
		t.Errorf("ExtractRegion failed: %v", err)
	}
	if !bytes.Equal(result.Letters, original.Letters) {
		t.Errorf("Expected original for empty regions")
	}
	
	// Test invalid region format
	_, err = ExtractRegion(original, "invalid")
	if err == nil {
		t.Error("Expected error for invalid region format")
	}
	
	// Test invalid indices
	_, err = ExtractRegion(original, "10:15")
	if err == nil {
		t.Error("Expected error for out of range indices")
	}
	
	_, err = ExtractRegion(original, "5:3")
	if err == nil {
		t.Error("Expected error for invalid index order")
	}
}

func TestNewReader(t *testing.T) {
	input := strings.NewReader("@test\nACGT\n+test\n!!!!\n")
	reader := NewReader(input)
	
	if reader == nil {
		t.Error("NewReader returned nil")
	}
}

func TestNewWriter(t *testing.T) {
	var buf bytes.Buffer
	writer := NewWriter(&buf)
	
	if writer == nil {
		t.Error("NewWriter returned nil")
	}
}

func TestNewScanner(t *testing.T) {
	input := strings.NewReader("@test\nACGT\n+test\n!!!!\n")
	reader := NewReader(input)
	scanner := NewScanner(reader)
	
	if scanner == nil {
		t.Error("NewScanner returned nil")
	}
}

func TestScanner(t *testing.T) {
	input := strings.NewReader("@test\nACGT\n+test\n!!!!\n")
	reader := NewReader(input)
	scanner := NewScanner(reader)
	
	// Test Next and Seq
	if !scanner.Next() {
		t.Error("Scanner.Next() returned false")
	}
	
	seq := scanner.Seq()
	if !bytes.Equal(seq.ID1, []byte("@test")) {
		t.Errorf("Expected ID1 @test, got %s", seq.ID1)
	}
	if !bytes.Equal(seq.Letters, []byte("ACGT")) {
		t.Errorf("Expected Letters ACGT, got %s", seq.Letters)
	}
	if !bytes.Equal(seq.ID2, []byte("+test")) {
		t.Errorf("Expected ID2 +test, got %s", seq.ID2)
	}
	if !bytes.Equal(seq.Quality, []byte("!!!!")) {
		t.Errorf("Expected Quality !!!!, got %s", seq.Quality)
	}
	
	// Test end of file
	if scanner.Next() {
		t.Error("Scanner.Next() returned true at end of file")
	}
	
	// Test error handling
	if err := scanner.Error(); err != nil && err != io.EOF {
		t.Errorf("Unexpected error: %v", err)
	}
}

func TestWriter(t *testing.T) {
	var buf bytes.Buffer
	writer := NewWriter(&buf)
	
	seq := Sequence{
		ID1:     []byte("@test"),
		Letters: []byte("ACGT"),
		ID2:     []byte("+test"),
		Quality: []byte("!!!!"),
	}
	
	n, err := writer.Write(seq)
	if err != nil {
		t.Errorf("Write failed: %v", err)
	}
	
	expected := "@test\nACGT\n+test\n!!!!\n"
	if buf.String() != expected {
		t.Errorf("Expected %q, got %q", expected, buf.String())
	}
	
	if n != len(expected) {
		t.Errorf("Expected %d bytes written, got %d", len(expected), n)
	}
}

func TestMaybeID1(t *testing.T) {
	if !maybeID1([]byte("@test")) {
		t.Error("maybeID1 should return true for line starting with @")
	}
	
	if maybeID1([]byte("test")) {
		t.Error("maybeID1 should return false for line not starting with @")
	}
	
	if maybeID1([]byte("")) {
		t.Error("maybeID1 should return false for empty line")
	}
}

func TestMaybeID2(t *testing.T) {
	if !maybeID2([]byte("+test")) {
		t.Error("maybeID2 should return true for line starting with +")
	}
	
	if maybeID2([]byte("test")) {
		t.Error("maybeID2 should return false for line not starting with +")
	}
	
	if maybeID2([]byte("")) {
		t.Error("maybeID2 should return false for empty line")
	}
}
