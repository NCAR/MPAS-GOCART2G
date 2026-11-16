#!/bin/ksh
# This script will submit a specified MPAS-GOCART forecasts

# CHANGE SECTION START - for adoption of this script by other users

n=0        # flag to determine if the job submissions will depend on previous jobs
year=2024  #$1    # First Input Argument
month=10   #$2   # Second Input Argument
day=15     # Third Input Argument
create_static=FALSE
run_nudging=FALSE
grid_name=x1.163842

proj=P19010000    # project key for Derecho submission

# Set the path of the script to its own location and display it
path_script=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
echo "The script is running from: $path_script"

path_mpas_src=/glade/work/$USER/MPAS_GOCART2G/gocartMPAS  # MPAS source code location
path_fcst=${path_script}/RUN_001							# Location where we want to store MPAS forecasts
path_gfs=/glade/campaign/acom/MUSICA/MPAS/input/Example-1/GFS			        # Location of GFS data for MPAS initialization
path_merra=/glade/campaign/acom/MUSICA/MPAS/input/Example-1/MERRA2   			# Location of MERRA-2 data for MPAS initialization
path_emis=/glade/campaign/acom/MUSICA/MPAS/input/Example-1/EMIS      			# Location of all emission files
path_back=/glade/campaign/acom/MUSICA/MPAS/input/Example-1/BACKGROUND			# Location of all background concentration files
path_grid=/glade/campaign/acom/MUSICA/MPAS/input/Example-1/GRID        			# Location of grid files for MPAS initialization
path_name=//glade/campaign/acom/MUSICA/MPAS/input/namelist_files			# Location of init_atmosphere namelist templates
path_stream=/glade/campaign/acom/MUSICA/MPAS/input/stream_files				# Location of init_atmosphere stream templates
fcst_len=00             	 # Forecast length in days (set to 00 for less than one day)
fcst_hr=03                       # Forecast length in hours (set to 00 for full days)
# CHANGE SECTION END

###
###
###
###

# Set MPAS fcst directory and define paths of namelist blueprints 
date_val=${year}${month}${day}
path_mpas_out=${path_fcst}/${date_val}
mkdir -p ${path_mpas_out}

# Build MPAS namelists on the fly: Calculate start and end dates and then execute the .sh files in fcst directory
start_year=$year
start_month=$month
start_day=$day
start_hour=00

echo "MPAS-GOCART-2G forecast start date is: ${start_year}-${start_month}-${start_day}_${start_hour}:00:00"

eventual_end_date=$(date +%C%y%m%d -d "$date_val +$fcst_len days")
end_year=`echo $eventual_end_date | cut -c1-4`
end_month=`echo $eventual_end_date | cut -c5-6`
end_day=`echo $eventual_end_date | cut -c7-8`
end_hour=$fcst_hr

echo "MPAS-GOCART-2G forecast end date is:   ${end_year}-${end_month}-${end_day}_${end_hour}:00:00"

render_stream_template() {
    stream_file=$1
    python3 - "$stream_file" "$grid_name" <<'PY'
import sys
from pathlib import Path
stream_path = Path(sys.argv[1])
grid_name = sys.argv[2]
text = stream_path.read_text()
text = text.replace('__GRID_NAME__', grid_name)
stream_path.write_text(text)
PY
}

prepare_init_atm_stage() {
    stage_name=$1
    stream_file=$2
    namelist_file=$3

    cd ${path_mpas_out}

    cp ${path_stream}/${stream_file} ./streams.${stage_name}.init_atmosphere
    cp ${path_name}/${namelist_file} ./${namelist_file}
    render_stream_template ./streams.${stage_name}.init_atmosphere

    . ./${namelist_file}
    cat << EOF > namelist.${stage_name}.init_atmosphere
$namelist_init_atmosphere
EOF
}

prepare_atm_stage() {
    stage_name=$1
    stream_file=$2
    namelist_file=$3

    cd ${path_mpas_out}

    cp ${path_stream}/stream_list* .
    cp ${path_stream}/${stream_file} ./streams.${stage_name}.atmosphere
    cp ${path_name}/${namelist_file} ./${namelist_file}
    render_stream_template ./streams.${stage_name}.atmosphere

    . ./${namelist_file}
    cat << EOF > namelist.${stage_name}.atmosphere
$namelist_atmosphere
EOF
}

link_common_files() {
    ln -sf ${path_grid}/${grid_name}.grid.nc .
    ln -sf ${path_grid}/${grid_name}.graph.info.part.* .
    ln -sf ${path_grid}/${grid_name}.static.nc .
    ln -sf ${path_grid}/${grid_name}.static_chems.nc .
    ln -sf ${path_back}/* .
    ln -sf ${path_emis}/*.nc .
    ln -sf ${path_gfs}/* .
    ln -sf ${path_merra}/* .
    ln -sf ${path_mpas_src}/atmosphere_model .
    ln -sf ${path_mpas_src}/init_atmosphere_model .
    ln -sf ${path_mpas_src}/*.TBL .
    ln -sf ${path_mpas_src}/*.DBL .
    ln -sf ${path_mpas_src}/CCN_ACTIVATE_DATA .
    ln -sf ${path_mpas_src}/RRTMG_LW_DATA .
    ln -sf ${path_mpas_src}/RRTMG_SW_DATA .
    ln -sf ${path_mpas_src}/src/core_atmosphere/physics/physics_wrf/files/COMPATIBILITY .
    ln -sf ${path_mpas_src}/src/core_atmosphere/physics/physics_wrf/files/VERSION .
    ln -sf /glade/campaign/acom/MUSICA/MPAS/input/optics_files/optics*.nc .
}

link_if_needed() {
    src=$1
    dst=$2
    if [ -n "${src}" ] && [ -n "${dst}" ] && [ "$(readlink -f "${src}" 2>/dev/null || echo "${src}")" != "$(readlink -f "${dst}" 2>/dev/null || echo "${dst}")" ]; then
        ln -sf "${src}" "${dst}"
    fi
}

write_init_atm_pbs() {
    run_label=$1
    walltime=$2
    output_file=$3
    namelist_name=$4
    stream_name=$5

    cat << EOF > ${output_file}
#!/bin/bash -l
#PBS -N ${run_label}
#PBS -A ${proj}
#PBS -q main
#PBS -l walltime=${walltime}
#PBS -m abe
#PBS -j oe
#PBS -l select=1:ncpus=128:mpiprocs=128:mem=235gb

mpiexec ./init_atmosphere_model -n ${namelist_name} -s ${stream_name}
EOF
}

write_atm_pbs() {
    run_label=$1
    walltime=$2
    output_file=$3
    namelist_name=$4
    stream_name=$5

    cat << EOF > ${output_file}
#!/bin/bash -l
#PBS -N ${run_label}
#PBS -A ${proj}
#PBS -q main
#PBS -l walltime=${walltime}
#PBS -m abe
#PBS -j oe
#PBS -l select=1:ncpus=128:mpiprocs=128:mem=235gb

mpiexec ./atmosphere_model -n ${namelist_name} -s ${stream_name}
EOF
}

cd ${path_mpas_out}
link_common_files

# Step 1: create static chemistry file if requested
if [ "${create_static}" = "TRUE" ]; then
    prepare_init_atm_stage "static" "streams.static.init_atmosphere" "namelist.static.init_atmosphere.sh"
    write_init_atm_pbs "init_static" "00:30:00" "${path_mpas_out}/run_init_static.pbs" "namelist.static.init_atmosphere" "streams.static.init_atmosphere"
    echo "Prepared static init_atmosphere run in: ${path_mpas_out}"
else
    echo "Skipping static init_atmosphere step because create_static=${create_static}"
fi

# Step 2: chemistry init_atmosphere stage, which reads the static file and produces init_chems_emissions.nc
prepare_init_atm_stage "chem" "streams.chem.init_atmosphere" "namelist.chem.init_atmosphere.sh"
if [ -f "${path_mpas_out}/${grid_name}.static_chems.nc" ]; then
    link_if_needed "${path_mpas_out}/${grid_name}.static_chems.nc" "./${grid_name}.static_chems.nc"
else
    echo "Warning: ${grid_name}.static_chems.nc not found in ${path_mpas_out}; chemistry stage may fail without the static file."
    echo "Create this file by changing create_static=TRUE and running run_init_static.pbs "
fi
write_init_atm_pbs "init_chem" "00:45:00" "${path_mpas_out}/run_init_chem.pbs" "namelist.chem.init_atmosphere" "streams.chem.init_atmosphere"
echo "Prepared chemistry init_atmosphere run in: ${path_mpas_out}"

# Step 3: nudging init_atmosphere stage, optional
if [ "${run_nudging}" = "TRUE" ]; then
    prepare_init_atm_stage "nudging" "streams.nudging.init_atmosphere" "namelist.nudging.init_atmosphere.sh"
    if [ -f "${path_mpas_out}/init_chems_emissions.nc" ]; then
        link_if_needed "${path_mpas_out}/init_chems_emissions.nc" "./init_chems_emissions.nc"
    else
        echo "Warning: init_chems_emissions.nc not found in ${path_mpas_out}; nudging stage may fail without the chemistry initialization."
        echo "Create this file by running run_init_chem.pbs."
    fi
    write_init_atm_pbs "init_nudging" "01:00:00" "${path_mpas_out}/run_init_nudging.pbs" "namelist.nudging.init_atmosphere" "streams.nudging.init_atmosphere"
    echo "Prepared nudging init_atmosphere run in: ${path_mpas_out}"

    prepare_atm_stage "nudging" "streams.nudging.atmosphere" "namelist.nudging.atmosphere.sh"
    if [ -f "${path_mpas_out}/${grid_name}.fdda.nc" ]; then
        link_if_needed "${path_mpas_out}/${grid_name}.fdda.nc" "./${grid_name}.fdda.nc"
    else
        echo "Warning: ${grid_name}.fdda.nc not found in ${path_mpas_out}; nudged atmosphere run may fail without FDDA input."
    fi
    if [ -f "${path_mpas_out}/init_chems_emissions.nc" ]; then
        link_if_needed "${path_mpas_out}/init_chems_emissions.nc" "./init_chems_emissions.nc"
    else
        echo "Warning: init_chems_emissions.nc not found in ${path_mpas_out}; atmosphere run may fail without the chemistry initialization."
        echo "Create this file by running run_init_chem.pbs."
    fi
    write_atm_pbs "atm_nudging" "01:00:00" "${path_mpas_out}/run_atm_nudging.pbs" "namelist.nudging.atmosphere" "streams.nudging.atmosphere"
    echo "Prepared nudged atmosphere run in: ${path_mpas_out}"
else
    echo "Skipping nudging init_atmosphere step because run_nudging=${run_nudging}"

    prepare_atm_stage "chem" "streams.chem.atmosphere" "namelist.chem.atmosphere.sh"
    if [ -f "${path_mpas_out}/init_chems_emissions.nc" ]; then
        link_if_needed "${path_mpas_out}/init_chems_emissions.nc" "./init_chems_emissions.nc"
    else
        echo "Warning: init_chems_emissions.nc not found in ${path_mpas_out}; chemistry atmosphere run may fail without the chemistry initialization."
        echo "Create this file by running run_init_chem.pbs."
    fi
    write_atm_pbs "atm_chem" "01:00:00" "${path_mpas_out}/run_atm_chem.pbs" "namelist.chem.atmosphere" "streams.chem.atmosphere"
    echo "Prepared chemistry atmosphere run in: ${path_mpas_out}"
fi

# Define jobs in strict execution order
job_scripts=(
    "run_init_static.pbs"
    "run_init_chem.pbs"
    "run_init_nudging.pbs"
    "run_atm_chem.pbs"
    "run_atm_nudging.pbs"
)

PREV_JOBID=""
FIRST_JOBID=""

for script in "${job_scripts[@]}"; do
    if [[ -f "$script" ]]; then
        if [[ -z "$PREV_JOBID" ]]; then
            # First existing job in chain: submit on hold (-h)
            echo "Submitting root job (held): $script"
            CURRENT_JOBID=$(qsub -h "$script") || exit 1
            FIRST_JOBID="$CURRENT_JOBID"
        else
            # Subsequent jobs: depend on the immediately preceding existing job
            echo "Submitting dependent job: $script (afterok:$PREV_JOBID)"
            CURRENT_JOBID=$(qsub -W depend=afterok:"$PREV_JOBID" "$script") || exit 1
        fi
        
        # Advance dependency pointer
        PREV_JOBID="$CURRENT_JOBID"
    else
        echo "Skipping $script (file not found)"
    fi
done

# Release the first job if at least one script was found and submitted
if [[ -n "$FIRST_JOBID" ]]; then
    echo "Releasing chain starting at $FIRST_JOBID"
    qrls "$FIRST_JOBID"
else
    echo "No matching .pbs files found. Nothing submitted."
fi

# Return to script directory
cd ${path_script}

# End of script 
exit 0
