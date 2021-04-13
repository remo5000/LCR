# Advogato
wget http://konect.cc/files/download.tsv.advogato.tar.bz2
bzip2 -df download.tsv.advogato.tar.bz2
tar -xf download.tsv.advogato.tar
cp advogato/out.advogato ./
python3 advogato.py
mv advogato.edge ../

rm download.tsv.advogato.tar
rm -rf advogato
rm out.advogato

# Epinions
wget http://konect.cc/files/download.tsv.soc-Epinions1.tar.bz2
bzip2 -df download.tsv.soc-Epinions1.tar.bz2
tar -xf download.tsv.soc-Epinions1.tar
cat soc-Epinions1/out.soc-Epinions1 | tail -n +3 > epinions.edge
python3 generateSyntheticLabels.py epinions.edge
mv epinions.edge ../

rm download.tsv.soc-Epinions1.tar
rm -rf soc-Epinions1

# Notre Dame
wget http://konect.cc/files/download.tsv.web-NotreDame.tar.bz2
bzip2 -df download.tsv.web-NotreDame.tar.bz2
tar -xf download.tsv.web-NotreDame.tar
cat web-NotreDame/out.web-NotreDame | tail -n +3 > web-NotreDame.edge
cat web-NotreDame.edge | ./generateSyntheticLabels.py > ../web-NotreDame.edge
mv web-NotreDame.edge ../

rm download.tsv.web-NotreDame.tar
rm -rf web-NotreDame/

# zhishi
wget http://konect.cc/files/download.tsv.zhishi-all.tar.bz2
bzip2 -df download.tsv.zhishi-all.tar.bz2
tar -xf download.tsv.zhishi-all.tar
cat zhishi-all/out.zhishi-all | tail -n +2 > zhishi-all.edge
python3 generateSyntheticLabels.py zhishi-all.edge
mv zhishi-all.edge ../

rm download.tsv.zhishi-all.tar
rm -rf zhishi-all

# wikilinks
wget http://konect.cc/files/download.tsv.wikipedia_link_fr.tar.bz2
bzip2 -df download.tsv.wikipedia_link_fr.tar.bz2
tar -xf download.tsv.wikipedia_link_fr.tar
cat wikipedia_link_fr/out.wikipedia_link_fr | tail -n +2 > wikipedia_link_fr.edge
python3 generateSyntheticLabels.py wikipedia_link_fr.edge
mv wikipedia_link_fr.edge ../

rm download.tsv.wikipedia_link_fr.tar
rm -rf wikipedia_link_fr

# webStanford
wget http://snap.stanford.edu/data/web-Stanford.txt.gz
gzip -d web-Stanford.txt.gz
tail -n +5 web-Stanford.txt > web-Stanford.edge
python3 generateSyntheticLabels.py web-Stanford.edge
mv web-Stanford.edge ../

rm web-Stanford.txt

# webGoogle
wget http://snap.stanford.edu/data/web-Google.txt.gz
gzip -d web-Google.txt.gz
tail -n +5 web-Google.txt > web-Google.edge
python3 generateSyntheticLabels.py web-Google.edge
mv web-Google.edge ../

rm web-Google.txt

#webBerkstan
wget http://snap.stanford.edu/data/web-BerkStan.txt.gz
gzip -d web-BerkStan.txt.gz
tail -n +5 web-BerkStan.txt > web-BerkStan.edge
python3 generateSyntheticLabels.py web-BerkStan.edge
mv web-BerkStan.edge ../

rm web-BerkStan.txt

#socPoekec
wget http://snap.stanford.edu/data/soc-pokec-relationships.txt.gz
gzip -d soc-pokec-relationships.txt.gz
mv soc-pokec-relationships.txt soc-pokec-relationships.edge
python3 generateSyntheticLabels.py soc-pokec-relationships.edge
mv soc-pokec-relationships.edge ../

# robots
wget http://tinyurl.com/gnexfoy -o robots_net-graph-2014-07-07.dot
# See robots.py on what to run in order to prepare the file for parsing.
# Then, run 
# python3 robots.py
mv robots.edge ../

# StringHS
wget http://version10.string-db.org/download/protein.actions.v10/9606.protein.actions.v10.txt.gz
gzip -d 9606.protein.actions.v10.txt.gz
python3 stringdb.py 9606.protein.actions.v10.txt stringHs.edge
mv stringHs.edge ../

rm 9606.protein.actions.v10.txt

# StringFC
wget http://version10.string-db.org/download/protein.actions.v10/9685.protein.actions.v10.txt.gz
gzip -d 9685.protein.actions.v10.txt.gz
python3 stringdb.py 9685.protein.actions.v10.txt stringFc.edge
mv stringFc.edge ../

rm 9685.protein.actions.v10.txt
