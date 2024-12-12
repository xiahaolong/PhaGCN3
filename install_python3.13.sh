cd Python-3.13.0
sudo ./configure --enable-optimizations --disable-gil --with-ssl
sudo make
sudo make altinstall
python3.13t -m pip install joblib
python3.13t -m pip install networkx
python3.13t -m pip install numpy
python3.13t -m pip install pandas

