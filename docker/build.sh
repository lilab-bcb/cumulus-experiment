docker build -t cumulus-experiment .
docker tag cumulus-experiment cumulusprod/cumulus-experiment
docker push cumulusprod/cumulus-experiment