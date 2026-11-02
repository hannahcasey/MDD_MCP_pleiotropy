## Start in empty project folder
mkdir resources
cd resources 

git clone https://github.com/comorment/mixer.git

mkdir oras
cd oras
wget https://github.com/oras-project/oras/releases/download/v1.3.0/oras_1.3.0_linux_amd64.tar.gz
tar -xzf oras_1.3.0_linux_amd64.tar.gz
chmod +x oras

cd ../mixer
../oras/oras pull ghcr.io/precimed/gsa-mixer_sif:2.2.1 && mv gsa-mixer.sif mixer.sif

mkdir 1000G
mkdir sumstats

cd ../../
mkdir output
mkdir bin