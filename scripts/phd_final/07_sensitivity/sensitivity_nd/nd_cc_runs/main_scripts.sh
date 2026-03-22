
# sequence 0500, 0600, ..., 1500 
vals=(
    0500
    0600
    0700
    0800
    0900
    1000
    1100
    1200
    1300
    1400
    1500
)

for val in "${vals[@]}"
do
    sbatch --export=STERILE=$val --job-name=$val run_nd.sh
done
