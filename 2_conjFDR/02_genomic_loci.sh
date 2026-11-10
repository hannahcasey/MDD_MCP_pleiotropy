cd python_convert

python3 -m venv venv
source venv/bin/activate
# pip install -U pip pandas numpy scipy

# wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/ref_1kG_phase3_EUR.tar.gz
# wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/9279485_ref.bim
# wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/2558411_ref.bim

# tar -xzf ref_1kG_phase3_EUR.tar.gz

PLINK=../plink_mac_20250615/plink
REFCHR=../python_convert/chr@

python sumstats.py clump \
  --clump-field FDR \
  --plink "$PLINK" \
  --sumstats conjFDR_MCP_MDD_for_clump.csv \
  --bfile-chr "$REFCHR" \
  --clump-p1 0.05 \
  --out ../../output/conjFDR/conjFDR_MCP_MDD.clump \
  --force