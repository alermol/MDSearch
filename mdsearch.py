#!/usr/bin/env python3
"""
MDSearch - Minimum Discriminatory SNPs set Search

A Python tool for identifying minimal sets of Single Nucleotide Polymorphisms (SNPs)
that can discriminate between samples in a VCF file based on a specified minimum Hamming distance.

This is the main entry point for the modular version of MDSearch.
"""

from src.cli import main

if __name__ == "__main__":
    main()
