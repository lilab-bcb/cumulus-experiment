docker build -t cumulus-experiment:19.11 .
docker tag cumulus-experiment:19.11 cumulusprod/cumulus-experiment:19.11
docker push cumulusprod/cumulus-experiment:19.11