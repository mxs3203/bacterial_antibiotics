#!/bin/bash

 # ./ARG_model.sh --convert Data-Fasta/Processed/ --load Data-Csv/ --test 10000 --model knn 1000 10 15 10
 # ./ARG_model.sh --convert Data-Fasta/Processed/

ITERATIONS=50
METHOD="GLM"
N_FEATURES=10
TUNE_LENGTH=15
CV_REPEATS=5
CSVFILE="test_predict.csv"

while test $# -gt 0; do   # check arguments one by one until there is none
  case "$1" in    # Perform actions for different cases
    -a|--all)
    shift         # remove first argument ($2 now is $1 if any) (remove flag (-a|--all))
    if test $# -gt 0; then    # get nº arguments after -l option (should be one -> CSVDIRECTORY)
      echo ""
      echo "Running all analysis with default values:"
      echo ""
      export CSVDIRECTORY=$1
      $0 -l $CSVDIRECTORY -t $ITERATIONS -m $METHOD $ITERATIONS $N_FEATURES $TUNE_LENGTH $CV_REPEATS 
      $0 -p $METHOD $CSVFILE 
      exit 0 
    else
      echo "No directory specified for loading CSV data"
      exit 1
    fi
    ;;

    -c|--convert)
    shift
    if test $# -gt 0; then    # get nº arguments after -c option (should be one -> CSVDIRECTORY)
      echo "Converting Fasta files in $1 directory to csv files"
      echo ""
      export FASTADIRECTORY=$1
      for filename in $FASTADIRECTORY/*; do
          python3 make_csv_from_fasta.py "$filename" "Data-Csv/$(basename "$filename" .fasta).csv"
      done
      shift
      # exit 0
    else
      echo "No directory specified for loading CSV data"
      exit 1
    fi
    ;;

    -s|--standardize)
    shift
    if test $# -gt 0; then    # get nº arguments after -c option (should be one -> CSVDIRECTORY)
      echo "Standardizing Fasta files in $1 directory and saving in Data-Fasta/Processed/"
      echo ""
      export FASTADIRECTORY=$1
      for filename in $FASTADIRECTORY/*; do
          python3 change_standard_of_fasta.py "$filename" "Data-Fasta/Processed/standardized$(basename "$filename" .fasta).fasta"
      done
      shift
      # exit 0
    else
      echo "No directory specified for converting fasta data"
      exit 1
    fi
    ;;

    -l|--load)
    shift         # remove first argument ($2 now is $1 if any) (rm "-l")
    if test $# -gt 0; then    # get nº arguments after -l option (should be one -> CSVDIRECTORY)
      echo "Load CSV Data into .RData file"
      echo ""
      export CSVDIRECTORY=$1
      Rscript LoadData.R $CSVDIRECTORY
      shift
      # exit 0
    else
      echo "No directory specified for loading CSV data"
      exit 1
    fi
    ;;
    
    -p|--predict)
    shift         # remove first argument ($2 now is $1 if any) (rm "-l")
    if test $# -gt 0; then   
      echo "Predicting antibiotic resistance genes"
      echo ""
      export CSVFILE=$2
      export METHOD=$1
      Rscript output_predictions.R $METHOD $CSVFILE
      shift 2
      # exit 0
    else
      echo "No File or method specified"
      exit 1
    fi
    ;;

    -t|--test)
    shift   # Remove flag (-t|--test)
    if test $# -gt 0; then    # get nº arguments after -t option (should be one -> ITERATIONS)
      echo "Find significantly different features between resistant and non-resistant genes:"
      echo ""
      export ITERATIONS=$1
      Rscript TestFeatures.R $ITERATIONS
      # exit 0
    else
      echo "Number of iterations not specified. Default nº of iterations = $ITERATIONS"
      Rscript TestFeatures.R $ITERATIONS
      # exit 0
    fi
    shift
    ;;

    -m|--model)
    shift
    if test $# -eq 5; then    # if nº arguments after -m option > 0
      export METHOD=$1
      export ITERATIONS=$2
      export N_FEATURES=$3
      export TUNE_LENGTH=$4
      export CV_REPEATS=$5
      shift 5
      Rscript BuildAModelAndPredict.R $METHOD $ITERATIONS $N_FEATURES $TUNE_LENGTH $CV_REPEATS
      
    else
      if test $# -eq 0; then
        echo "Fitting and testing the model with default values"
        Rscript BuildAModelAndPredict.R $METHOD $ITERATIONS $N_FEATURES $TUNE_LENGTH $CV_REPEATS
        
      fi
      echo "Incorrect number of arguments specified. Five or zero arguments expected"
      echo ""
      echo "Usage: $0 -m [METHOD] [ITERATIONS] [N_FEATURES] [TUNE_LENGTH] [CV_REPEATS]"
      echo "Usage (defaults): $0 -m"
      exit 1
    fi
    ;;

    -h|--help)
    echo "$0:    ANTIBIOTIC RESISTANT GENES MODEL"
    echo ""
    echo "You should follow the next pipeline for the analysis:"
    echo "    1) Load the data (option '-l' or '--load')"
    echo "    2) Find significantly different features (option '-t' or '--test')"
    echo "    3) Create model and test predictions (option '-m' or '--model')"
    echo ""
    echo "OR use '$0 -a [CSVDIRECTORY]' command to run whole pipeline with default parameters"
    echo ""
    echo " "
    echo "$0 [options] [arguments]"
    echo " "
    echo "Options:"
    echo ""
    echo "  -a, --all                 run all analysis with default values"
    echo "  -c, --convert             convert fasta files to csv files"
    echo "  -h, --help                display this help"
    echo "  -l, --load                load csv files to .RData."
    echo "                              CSVDIRECTORY should be input as argument"
    echo "  -m, --model               create model and give test results back."
    echo "                              Several arguments can be added."
    echo "                              ITERATIONS should be input as argument"
    # echo "  -o, --output-dir=DIR      specify a directory to store output in"
    echo "  -p, --predict                predict antibiotic resistance genes from csv file"
    echo "  -s, --standardize         standardize fasta files for analysis"
    echo "  -t, --test                test for significantly different features"
    echo "                              between resistant and non-resistant genes."
    echo "                              ITERATIONS should be input as argument"
    echo ""
    echo "Arguments:"
    echo ""
    echo "  -c | -s option:"
    echo "      FASTADIRECTORY"
    echo ""
    echo "  -l | -a options:"
    echo "      CSVDIRECTORY            directory with the CSV files"
    echo ""
    echo "  -t option:"
    echo "      ITERATIONS              number of iterations used in resampling"
    echo ""
    echo "  -m option:"
    echo "      METHOD                 statistical method used to fit the model (KNN/GLM)"
    echo "      ITERATIONS             number of iterations used to fit and test the model"
    echo "      N_FEATURES             number of (most significant) features used to fit the model"
    echo "      TUNE_LENGTH            KNN tune length"
    echo "      CV_REPEATS             number of repeats used in cross validation when method = KNN"
    echo "  *arguments for -m should be input in the next order:"
    echo "    -m [METHOD] [ITERATIONS] [N_FEATURES] [TUNE_LENGTH] [CV_REPEATS]"
    echo "" 
    echo "  -p option:"
    echo "      METHOD                 statistical method used to fit the model (KNN/GLM)"
    echo "      CSVFILE                path to CSV file"
    exit 0
    ;;

    *)    # other arguments
    echo ""
    echo "Yo no comprendo :("
    echo ""
    echo "Please use '$0 --help' for usage info"
    exit 1
    ;;
  esac

done
