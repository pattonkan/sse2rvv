package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

func main() {
	sourceFileName := "../sse2rvv.h"
	modifiedFileName := "../modified_sse2rvv.h"
	cleanUpFile(sourceFileName, modifiedFileName)
	os.Remove(sourceFileName)
	os.Rename(modifiedFileName, sourceFileName)
}

func cleanUpFile(sourceFileName string, modifiedFileName string) {
	file, err := os.Open(sourceFileName)
	if err != nil {
		fmt.Printf("Error opening file: %v\n", err)
		return
	}
	defer file.Close()

	modifiedFile, err := os.Create(modifiedFileName)
	if err != nil {
		fmt.Printf("Error creating modified file: %v\n", err)
		return
	}
	defer modifiedFile.Close()

	scanner := bufio.NewScanner(file)
	commentType := 0
	bracket := 0
	for scanner.Scan() {
		line := scanner.Text()
		if getCommentType(line) > 0 {
			commentType = getCommentType(line)
		}

		if commentType > 0 {
			if commentType == 1 {
				if strings.ContainsAny(line, "{") {
					bracket += 1
				}
				if strings.ContainsAny(line, "}") {
					bracket -= 1
				}

				trimmed := strings.TrimSpace(line)
				if strings.HasPrefix(trimmed, "FORCE_INLINE") && strings.HasSuffix(trimmed, ";") || strings.HasSuffix(trimmed, ",") && bracket == 0 {
					modifiedFile.WriteString("// " + line + "\n")
				} else if strings.HasPrefix(trimmed, "FORCE_INLINE") || strings.HasSuffix(trimmed, ")") && bracket == 0 {
					modifiedFile.WriteString("// " + line + " ")
				} else if trimmed == "{" {
					modifiedFile.WriteString(trimmed)
				} else if trimmed == "}" && bracket == 0 {
					modifiedFile.WriteString(trimmed + "\n")
				}

				if bracket > 0 {
					continue
				}

				if strings.HasSuffix(trimmed, "}") || strings.HasSuffix(trimmed, ";") {
					commentType = 0
				}
			} else if commentType == 2 {
				trimmed := strings.TrimSpace(line)
				if strings.HasPrefix(trimmed, "#define") {
					s := strings.ReplaceAll(line, "\\", "")
					modifiedFile.WriteString("// " + strings.TrimSpace(s) + "\n")
				}

				if !strings.HasSuffix(trimmed, "\\") {
					commentType = 0
				}
			}
		} else {
			modifiedFile.WriteString(line + "\n")
		}
	}

	if err := scanner.Err(); err != nil {
		fmt.Printf("Error reading file: %v\n", err)
		return
	}
}

func getCommentType(line string) int {
	trimmed := strings.TrimSpace(line)
	if strings.HasPrefix(trimmed, "FORCE_INLINE") && strings.Contains(trimmed, "_mm_") {
		return 1
	} else if strings.HasPrefix(trimmed, "#define") && strings.Contains(trimmed, "_mm_") {
		return 2
	} else {
		return 0
	}
}
