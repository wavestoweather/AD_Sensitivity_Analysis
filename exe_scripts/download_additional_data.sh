#!/bin/bash
cd ..
if [ ! -d "data" ]
then
    mkdir -p "data"
fi
cd "data"
echo "Information about the data to be downloaded:"
curl https://irods-web.zdv.uni-mainz.de/irods-rest/rest/dataObject/zdv/project/m2_jgu-w2w/w2w-z2/vladiana_complete.tar.gz/metadata?ticket=WAHrBhMLpBn5ZVT
wget https://irods-web.zdv.uni-mainz.de/irods-rest/rest/fileContents/zdv/project/m2_jgu-w2w/w2w-z2/vladiana_complete.tar.gz?ticket=WAHrBhMLpBn5ZVT -O vladiana_complete.tar.gz
echo "Downloaded data. Now extracting."
tar -xvf vladiana_complete.tar.gz
rm vladiana_complete.tar.gz
echo "Data has been downloaded and extracted to data/vladiana_complete/"