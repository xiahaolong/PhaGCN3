使用方法
git clone https://github.com/xiahaolong/VirusGCN.git
cd VirusGCN
conda env create -f VirusGCN.yaml -n virusgcn
安装python自由线程版本python3.13t
构建安装环境和下载python3.13
sudo apt-get upgrade
sudo apt-get install build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev
sudo wget https://www.python.org/ftp/python/3.13.0/Python-3.13.0.tgz
sudo tar xzf Python-3.13.0.tgz
脚本安装python3.13，以及需要的包
./install_python3.13.sh
数据库解压
cd database
tar -zxvf ALL_protein.tar.gz
cd ..
运行程序
conda activate virusgcn
python run_Speed_up.py --contigs contigs.fa --len 8000 --outpath result2