"""Tests for SNP Profile Writer module."""

import pytest
from pathlib import Path
from unittest.mock import Mock, MagicMock

from src.io.snp_profile_writer import SNPProfileWriter
from src.core.vcf_parser import VCFData, VCFHeaders, SNPData


class TestSNPProfileWriter:
    """Test cases for SNPProfileWriter class."""

    def test_init_with_logger(self):
        """Test initialization with custom logger."""
        mock_logger = Mock()
        writer = SNPProfileWriter(mock_logger)
        assert writer.logger == mock_logger

    def test_init_without_logger(self):
        """Test initialization without custom logger."""
        writer = SNPProfileWriter()
        assert writer.logger is not None

    def test_write_snp_profiles_creates_directory(self, tmp_path):
        """Test that write_snp_profiles creates the snp_profiles directory."""
        writer = SNPProfileWriter()
        
        # Create mock VCF data
        mock_vcf_data = Mock(spec=VCFData)
        mock_vcf_data.headers = Mock()
        mock_vcf_data.headers.samples = ["sample1", "sample2"]
        mock_vcf_data.snp_genotypes = {}
        
        snp_sets = [["rs1", "rs2"]]
        
        profiles_dir = writer.write_snp_profiles(tmp_path, snp_sets, mock_vcf_data)
        
        assert profiles_dir.exists()
        assert profiles_dir.name == "snp_profiles"

    def test_write_snp_profiles_generates_files(self, tmp_path):
        """Test that write_snp_profiles generates profile files for each SNP set."""
        writer = SNPProfileWriter()
        
        # Create mock VCF data with SNP information
        mock_vcf_data = Mock(spec=VCFData)
        mock_vcf_data.headers = Mock()
        mock_vcf_data.headers.samples = ["sample1", "sample2"]
        mock_vcf_data.snp_genotypes = {
            "rs1": Mock(spec=SNPData),
            "rs2": Mock(spec=SNPData)
        }
        
                # Configure mock SNP data
        mock_vcf_data.snp_genotypes["rs1"].chromosome = "chr1"
        mock_vcf_data.snp_genotypes["rs1"].position = 1000
        mock_vcf_data.snp_genotypes["rs1"].maf = 0.3
        mock_vcf_data.snp_genotypes["rs1"].genotypes = [0, 1]
        mock_vcf_data.snp_genotypes["rs1"].reference_allele = "A"
        mock_vcf_data.snp_genotypes["rs1"].alternate_allele = "T"
        mock_vcf_data.snp_genotypes["rs1"].sample_fields = ["0/0", "0/1"]
    
        mock_vcf_data.snp_genotypes["rs2"].chromosome = "chr2"
        mock_vcf_data.snp_genotypes["rs2"].position = 2000
        mock_vcf_data.snp_genotypes["rs2"].maf = 0.4
        mock_vcf_data.snp_genotypes["rs2"].genotypes = [1, 0]
        mock_vcf_data.snp_genotypes["rs2"].reference_allele = "C"
        mock_vcf_data.snp_genotypes["rs2"].alternate_allele = "G"
        mock_vcf_data.snp_genotypes["rs2"].sample_fields = ["1/1", "0/0"]
        
        mock_vcf_data.snp_entropy_cache = {"rs1": 0.5, "rs2": 0.6}
        
        snp_sets = [["rs1", "rs2"]]
        
        profiles_dir = writer.write_snp_profiles(tmp_path, snp_sets, mock_vcf_data)
        
        # Check that profile file was created
        profile_file = profiles_dir / "snp_set_1_profile.txt"
        assert profile_file.exists()
        
        # Check file content
        content = profile_file.read_text()
        assert "sample1:" in content
        assert "sample2:" in content
        assert "chr1:1000" in content
        assert "chr2:2000" in content

    def test_write_snp_profiles_with_multiple_sets(self, tmp_path):
        """Test that write_snp_profiles handles multiple SNP sets."""
        writer = SNPProfileWriter()
        
        # Create mock VCF data
        mock_vcf_data = Mock(spec=VCFData)
        mock_vcf_data.headers = Mock()
        mock_vcf_data.headers.samples = ["sample1"]
        mock_vcf_data.snp_genotypes = {
            "rs1": Mock(spec=SNPData),
            "rs2": Mock(spec=SNPData)
        }
        
        # Configure mock SNP data
        for snp_id in ["rs1", "rs2"]:
            mock_vcf_data.snp_genotypes[snp_id].chromosome = "chr1"
            mock_vcf_data.snp_genotypes[snp_id].position = 1000
            mock_vcf_data.snp_genotypes[snp_id].maf = 0.3
            mock_vcf_data.snp_genotypes[snp_id].genotypes = [0]
            mock_vcf_data.snp_genotypes[snp_id].reference_allele = "A"
            mock_vcf_data.snp_genotypes[snp_id].alternate_allele = "T"
            mock_vcf_data.snp_genotypes[snp_id].sample_fields = ["0/0"]
        
        mock_vcf_data.snp_entropy_cache = {"rs1": 0.5, "rs2": 0.6}
        
        snp_sets = [["rs1"], ["rs2"]]
        
        profiles_dir = writer.write_snp_profiles(tmp_path, snp_sets, mock_vcf_data)
        
        # Check that both profile files were created
        assert (profiles_dir / "snp_set_1_profile.txt").exists()
        assert (profiles_dir / "snp_set_2_profile.txt").exists()

    def test_load_recode_scheme_valid_file(self, tmp_path):
        """Test loading recoding scheme from valid TSV file."""
        writer = SNPProfileWriter()
        
        # Create test recode file
        recode_file = tmp_path / "recode.tsv"
        recode_file.write_text("old1\tnew1\nold2\tnew2\n")
        
        scheme = writer._load_recode_scheme(recode_file)
        
        assert scheme["old1"] == "new1"
        assert scheme["old2"] == "new2"

    def test_load_recode_scheme_invalid_file(self, tmp_path):
        """Test loading recoding scheme from invalid file."""
        writer = SNPProfileWriter()
        
        # Create invalid recode file
        recode_file = tmp_path / "invalid.tsv"
        recode_file.write_text("invalid_format\n")
        
        scheme = writer._load_recode_scheme(recode_file)
        
        assert scheme == {}

    def test_recode_samples(self, tmp_path):
        """Test sample name recoding functionality."""
        writer = SNPProfileWriter()
        
        # Create test recode file
        recode_file = tmp_path / "samples.tsv"
        recode_file.write_text("sample1\tSAMPLE_001\nsample2\tSAMPLE_002\n")
        
        original_names = ["sample1", "sample2", "sample3"]
        recoded_names = writer._recode_samples(original_names, recode_file)
        
        assert recoded_names == ["SAMPLE_001", "SAMPLE_002", "sample3"]
