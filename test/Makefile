test:
	python3 generate_vcf.py test_ploidy2.vcf --ploidy 2 --nsamples 1000
	python3 ../mdsearch.py test_ploidy2.vcf out_pl2.vcf -e 100 -t 10
	python3 generate_vcf.py test_ploidy4.vcf --ploidy 4 --nsamples 1000
	python3 ../mdsearch.py test_ploidy4.vcf out_pl4.vcf -e 100 -t 10 -pl 4
	python3 generate_vcf.py test_ploidy6.vcf --ploidy 6 --nsamples 1000
	python3 ../mdsearch.py test_ploidy6.vcf out_pl6.vcf -e 100 -t 10 -pl 6
	python3 generate_vcf.py test_ploidy2_phased.vcf --ploidy 2 --phased --nsamples 1000
	python3 ../mdsearch.py test_ploidy2_phased.vcf out_pl2_phased.vcf -e 100 -t 10
	python3 generate_vcf.py test_ploidy4_phased.vcf --ploidy 4 --phased --nsamples 1000
	python3 ../mdsearch.py test_ploidy4_phased.vcf out_pl4_phased.vcf -e 100 -t 10 -pl 4
	python3 generate_vcf.py test_ploidy6_phased.vcf --ploidy 6 --phased --nsamples 1000
	python3 ../mdsearch.py test_ploidy6_phased.vcf out_pl6_phased.vcf -e 100 -t 10 -pl 6

clean:
	rm -v test*.vcf out*.vcf
