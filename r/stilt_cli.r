#!/usr/bin/env Rscript
# STILT R CLI for single-shot batch proccessing
# For documentation, see https://uataq.github.io/stilt/
# Ben Fasoli
# 
# git clone https://github.com/uataq/stilt-tutorials /tmp/stilt-tutorials
# r/run_stilt_cli.r \
#     r_run_time=2015-12-10T00:00:00Z \
#     r_lati=40.5 \
#     r_long=-112.0 \
#     r_zagl=5 \
#     met_loc=/tmp/stilt-tutorials/01-wbb/met/ \
#     met_file_format=%Y%m%d.%H \
#     xmn=-112.3 \
#     xmx=-111.52 \
#     xres=0.01 \
#     ymn=40.39 \
#     ymx=40.95 \
#     yres=0.01

# Extract kv pairs for supplied arguments -------------------------------------
arg_strings <- commandArgs(trailingOnly = T)
args <- list()
for (arg in strsplit(arg_strings, '=', fixed = T)) {
    if (length(arg) == 1) {
        arg <- c(arg, 'NA')
    }
    args[arg[1]] <- paste(arg[2:length(arg)], collapse='=')
}

# Validate required arguments exist
req_args <- c('met_file_format', 'met_loc', 'r_run_time', 'r_lati', 'r_long',
              'r_zagl', 'xmn', 'xmx', 'xres', 'ymn', 'ymx', 'yres')
if (!all(req_args %in% names(args))) {
    stop(paste('Not all arguments supplied:', 
                paste(req_args, collapse=',')))
}

# Ensure script is executed from the correct place
for (i in 1:99) {
    if ('stilt_wd' %in% names(args)) break
    cwd <- dir()
    is_stilt <- all(c('exe', 'fortran', 'r', 'README.md', 'setup') %in% cwd)
    if (is_stilt) {
        args$stilt_wd <- getwd()
    } else {
        setwd('..')
    }
}
if (i == 99) {
    stop('Could not identify STILT project directory - specify stilt_wd manually')
}

# Set argument types
stilt_args <- list(
    capemin <- as.numeric(args$capemin),
    cmass <- as.numeric(args$cmass),
    conage = as.numeric(args$conage),
    cpack = as.numeric(args$cpack),
    dxf = as.numeric(args$dxf),
    dyf = as.numeric(args$dyf),
    dzf = as.numeric(args$dzf),
    efile = as.character(args$efile),
    emisshrs = as.numeric(args$emisshrs),
    frhmax = as.numeric(args$frhmax),
    frhs = as.numeric(args$frhs),
    frme = as.numeric(args$frme),
    frmr = as.numeric(args$frmr),
    frts = as.numeric(args$frts),
    frvs = as.numeric(args$frvs),
    hnf_plume = as.logical(args$hnf_plume),
    horcoruverr = as.numeric(args$horcoruverr),
    horcorzierr = as.numeric(args$horcorzierr),
    hscale = as.numeric(args$hscale),
    ichem = as.numeric(args$ichem),
    iconvect = as.numeric(args$iconvect),
    idsp = as.numeric(args$idsp),
    initd = as.numeric(args$initd),
    isot = as.numeric(args$isot),
    k10m = as.numeric(args$k10m),
    kagl = as.numeric(args$kagl),
    kbls = as.numeric(args$kbls),
    kblt = as.numeric(args$kblt),
    kdef = as.numeric(args$kdef),
    khinp = as.numeric(args$khinp),
    khmax = as.numeric(args$khmax),
    kmix0 = as.numeric(args$kmix0),
    kmixd = as.numeric(args$kmixd),
    kmsl = as.numeric(args$kmsl),
    kpuff = as.numeric(args$kpuff),
    krand = as.numeric(args$krand),
    krnd = as.numeric(args$krnd),
    kspl = as.numeric(args$kspl),
    kwet = as.numeric(args$kwet),
    kzmix = as.numeric(args$kzmix),
    lib.loc = as.character(args$lib.loc),
    maxdim = as.numeric(args$maxdim),
    maxpar = as.numeric(args$maxpar),
    met_file_format = as.character(args$met_file_format),
    met_loc = as.character(args$met_loc),
    mgmin = as.numeric(args$mgmin),
    n_hours = as.numeric(args$n_hours),
    n_met_min = as.numeric(args$n_met_min),
    ncycl = as.numeric(args$ncycl),
    ndump = as.numeric(args$ndump),
    ninit = as.numeric(args$ninit),
    nstr = as.numeric(args$nstr),
    nturb = as.numeric(args$nturb),
    numpar = as.numeric(args$numpar),
    nver = as.numeric(args$nver),
    outdt = as.numeric(args$outdt),
    outfrac = as.numeric(args$outfrac),
    output_wd = as.character(args$output_wd),
    p10f = as.numeric(args$p10f),
    pinbc = as.character(args$pinbc),
    pinpf = as.character(args$pinpf),
    poutf = as.character(args$poutf),
    projection = as.character(args$projection),
    qcycle = as.numeric(args$qcycle),
    r_run_time = as.POSIXct(args$r_run_time,
                            tz = 'UTC',
                            format = '%Y-%m-%dT%H:%M:%SZ'),
    r_lati = as.numeric(args$r_lati),
    r_long = as.numeric(args$r_long),
    r_zagl = as.numeric(args$r_zagl),
    random = as.numeric(args$random),
    rhb = as.numeric(args$rhb),
    rht = as.numeric(args$rht),
    rm_dat = as.logical(args$rm_dat),
    run_foot = as.logical(args$run_foot),
    run_trajec = as.logical(args$run_trajec),
    siguverr = as.numeric(args$siguverr),
    sigzierr = as.numeric(args$sigzierr),
    smooth_factor = as.numeric(args$smooth_factor),
    splitf = as.numeric(args$splitf),
    stilt_wd = as.character(args$stilt_wd),
    time_integrate = as.logical(args$time_integrate),
    timeout = as.numeric(args$timeout),
    tkerd = as.numeric(args$tkerd),
    tkern = as.numeric(args$tkern),
    tlfrac = as.numeric(args$tlfrac),
    tluverr = as.numeric(args$tluverr),
    tlzierr = as.numeric(args$tlzierr),
    tout = as.numeric(args$tout),
    tratio = as.numeric(args$tratio),
    tvmix = as.numeric(args$tvmix),
    varsiwant = as.character(args$varsiwant),
    veght = as.numeric(args$veght),
    vscale = as.numeric(args$vscale),
    vscaleu = as.numeric(args$vscaleu),
    vscales = as.numeric(args$vscales),
    w_option = as.numeric(args$w_option),
    wbbh = as.numeric(args$wbbh),
    wbwf = as.numeric(args$wbwf),
    wbwr = as.numeric(args$wbwr),
    wvert = as.logical(args$wvert),
    xmn = as.numeric(args$xmn),
    xmx = as.numeric(args$xmx),
    xres = as.numeric(args$xres),
    ymn = as.numeric(args$ymn),
    ymx = as.numeric(args$ymx),
    yres = as.numeric(args$yres),
    zicontroltf = as.numeric(args$zicontroltf),
    ziscale = as.numeric(args$ziscale),
    z_top = as.numeric(args$z_top),
    zcoruverr = as.numeric(args$zcoruverr)
)
stilt_args <- stilt_args[sapply(stilt_args, function(x) length(x) > 0)]

source(file.path(stilt_args$stilt_wd, 'r', 'src', 'simulation_step.r'))
res <- do.call(simulation_step, stilt_args)
q('no')
