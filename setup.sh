mkdir setup &&
cd setup &&
wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py &&
python3 install.py --deps --tool --ot --zk &&
cd .. &&
sudo apt install -y emacs iperf iftop clang

