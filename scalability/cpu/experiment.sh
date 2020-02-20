BLUE='\033[1;36m'
NC='\033[0m'

for cores in 1 3 7 14
do
    echo -e "${BLUE}Pegasus experiment using ${cores} cores:${NC}"
    python run_pegasus.py $cores
    echo -e "${BLUE}SCANPY experiment using ${cores} cores:${NC}"
    python run_scanpy.py $cores > "scanpy_${cores}_cores.log"
    echo -e "${BLUE}Seurat experiment using ${cores} cores:${NC}"
    Rscript run_seurat.R --args $cores
done