Rscript --vanilla 1_run_sparse_model.R original 5 1 12141 &> progress/sparse_idfOrig_K5_inf1_12141.txt &
Rscript --vanilla 1_run_sparse_model.R original 5 1 12142 &> progress/sparse_idfOrig_K5_inf1_12142.txt &

Rscript --vanilla 1_run_sparse_model.R original 5 3 12141 &> progress/sparse_idfOrig_K5_inf3_12141.txt &
Rscript --vanilla 1_run_sparse_model.R original 5 3 12142 &> progress/sparse_idfOrig_K5_inf3_12142.txt &

Rscript --vanilla 1_run_sparse_model.R original 7 1 12141 &> progress/sparse_idfOrig_K7_inf1_12141.txt &
Rscript --vanilla 1_run_sparse_model.R original 7 1 12142 &> progress/sparse_idfOrig_K7_inf1_12142.txt &

Rscript --vanilla 1_run_sparse_model.R original 7 3 12141 &> progress/sparse_idfOrig_K7_inf3_12141.txt &
Rscript --vanilla 1_run_sparse_model.R original 7 3 12142 &> progress/sparse_idfOrig_K7_inf3_12142.txt &

# Rscript --vanilla 1_run_sparse_model.R new1 5 1 12141 &> progress/sparse_idfNew1_K5_inf1_12141.txt &
# Rscript --vanilla 1_run_sparse_model.R new1 5 3 12141 &> progress/sparse_idfNew1_K5_inf3_12141.txt &
# 
# Rscript --vanilla 1_run_sparse_model.R new1 7 1 12141 &> progress/sparse_idfNew1_K7_inf1_12141.txt &
# Rscript --vanilla 1_run_sparse_model.R new1 7 3 12141 &> progress/sparse_idfNew1_K7_inf3_12141.txt &


Rscript --vanilla 1_run_base_model.R original 5 1 12141 &> progress/base_idfOrig_K5_inf1_12141.txt &
Rscript --vanilla 1_run_base_model.R original 5 1 12142 &> progress/base_idfOrig_K5_inf1_12142.txt &

Rscript --vanilla 1_run_base_model.R original 5 3 12141 &> progress/base_idfOrig_K5_inf3_12141.txt &
Rscript --vanilla 1_run_base_model.R original 5 3 12142 &> progress/base_idfOrig_K5_inf3_12142.txt &

Rscript --vanilla 1_run_base_model.R original 7 1 12141 &> progress/base_idfOrig_K7_inf1_12141.txt &
Rscript --vanilla 1_run_base_model.R original 7 1 12142 &> progress/base_idfOrig_K7_inf1_12142.txt &

Rscript --vanilla 1_run_base_model.R original 7 3 12141 &> progress/base_idfOrig_K7_inf3_12141.txt &
Rscript --vanilla 1_run_base_model.R original 7 3 12142 &> progress/base_idfOrig_K7_inf3_12142.txt &

# Rscript --vanilla 1_run_base_model.R new1 5 1 12141 &> progress/base_idfNew1_K5_inf1_12141.txt &
# Rscript --vanilla 1_run_base_model.R new1 5 3 12141 &> progress/base_idfNew1_K5_inf3_12141.txt &
# 
# Rscript --vanilla 1_run_base_model.R new1 7 1 12141 &> progress/base_idfNew1_K7_inf1_12141.txt &
# Rscript --vanilla 1_run_base_model.R new1 7 3 12141 &> progress/base_idfNew1_K7_inf3_12141.txt &

wait
