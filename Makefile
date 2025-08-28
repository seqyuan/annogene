.PHONY: help build test clean lint format push

# Default target
help:
	@echo "Available commands:"
	@echo "  build    - Build the package"
	@echo "  test     - Run tests"
	@echo "  clean    - Clean build artifacts"
	@echo "  lint     - Run linter"
	@echo "  format   - Format code"
	@echo "  push     - Git add, commit and push"

# Build the package
build:
	go build ./...

# Run tests
test:
	go test -v ./...

# Clean build artifacts
clean:
	go clean
	rm -f *.exe *.test

# Run linter
lint:
	golangci-lint run

# Format code
format:
	go fmt ./...
	gofmt -s -w .

# Git operations
push:
	git add -A
	git commit -m "Update package"
	git push
	go get -u github.com/seqyuan/annogene

# Install dependencies
deps:
	go mod tidy
	go mod download

# Check for security vulnerabilities
security:
	govulncheck ./...
