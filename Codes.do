// Codes of Study "..."
//
// 1. Define .shp file for Spatial Weight Matrix
*  it could be used https://mapshaper.org/ 
*  After you get all world country map, use the web tool above and extract
*  whatever countries you want. 
*  IMPORTANT!: When you load data to "mapshaper" load it as .zip file!
*  Select the countries according to the economic data panel id order
*  Extract it as .zip file and unzip it, copy .shp and .dbf files to working direction
*  Then read it from Stata and get CX & CY values.
*  Your economic data order of countries will not match with .shp file.
*  Change ID variable according to your economic data set country order.

// read bri.shp & bri.dbf file and get Stata data file bri_data.dta & bri_coord.dta
// before you do that, you need to import .dbf file and export, save it again using
// Stata. Because the "shp2dta" command does not work Stata 15 and above.
import dbase using bri.dbf, clear
export dbase bri.dbf, replace
// now you could use "shp2dta" command.
shp2dta using bri, database(bri) coordinates(bri_xy) genid(id) gencentroids(c)
use bri
gen st = id
rename x_c longitude
rename y_c latitude
save bri, replace
// check the map
spmap using bri_xy, id(id)
// make the weight matrix
spmat contiguity M using bri_xy, id(id) normalize(row)
spmat summarize M
// load panel economic data
use paneldata, clear
xtset id year
// then estimate Dynamic Panel Spatial Regressions
// specification process: first we need to decide fe or re in SDM.
// because SDM is most inclusive model of all spatial models
* LM  Standard Panel Model Estimation
xtreg lnpy lnc lnh lnex lnim lnfd, fe  // FE
estimates store fixed
*
xtreg lnpy lnc lnh lnex lnim lnfd, re  // RE
estimates store random
*
hausman fixed random
* SDM
xsmle lnpy lnc lnh lnex lnim lnfd dbri, wmat(M) model(sdm) fe type(ind) durbin(lnc lnh lnim lnex lnfd) nolog // FE
estimates store sdm_fe
*
xsmle lnpy lnc lnh lnex lnim lnfd dbri, wmat(M) model(sdm) re type(ind) durbin(lnc lnh lnim lnex lnfd) nolog // RE
estimates store sdm_re
*
hausman sdm_fe sdm_re, eq(1:1 2:2 3:3)
* SAR
xsmle lnpy lnc lnh lnex lnim lnfd dbri, wmat(M) model(sar) fe type(ind) nolog
* SAC
xsmle lnpy lnc lnh lnex lnim lnfd dbri, wmat(M) emat(M) model(sac) fe type(ind) nolog
* SEM
xsmle lnpy lnc lnh lnex lnim lnfd dbri, emat(M) model(sem) fe type(ind) nolog

// Mapping 
use bri, clear
// GDP 2004
spmap y2004 using "bri_xy.dta", id(id) clmethod(custom) clbreaks(0 2000 4000 6000 8000) fcolor(Blues) legend(size(*1.5) position (4)) legorder(hilo) legend(label(2 "0 -  $2000") label(3 "$2000 - $4000") label(4 "$4000 - $6000") label(5 "$6000 - $8000")) label(xcoord(longitude) ycoord(latitude) label(y2004))
 
// GDP 2020
spmap y2020 using "bri_xy.dta", id(id) clmethod(custom) clbreaks(0 2000 4000 6000 8000 10000 12000) fcolor(Blues) legend(size(*1.5) position (4)) legorder(hilo) legend(label(2 "0 -  $2000") label(3 "$2000 - $4000") label(4 "$4000 - $6000") label(5 "$6000 - $8000") label(6 "$8000 - $10000") label(7 "$10000 - $12000")) label(xcoord(longitude) ycoord(latitude) label(y2020))

// FDI 2004
gen F2004=fd2004/1000000

spmap F2004 using "bri_xy.dta", id(id) clmethod(custom) clbreaks(0 10 20 30 40 50) fcolor(Blues) legend(size(*1.5) position (4)) legorder(hilo) legend(label(2 "0 -  $10") label(3 "$10 - $20") label(4 "$20 - $30") label(5 "$30 - $40") label(6 "$40 - $50")) label(xcoord(longitude) ycoord(latitude) label(F2004))

// FDI 2020
gen F2020=fd2020/1000000

spmap F2020 using "bri_xy.dta", id(id) clmethod(custom) clbreaks(0 1000 2000 3000 4000 5000 6000) fcolor(Blues) legend(size(*1.5) position (4)) legorder(hilo) legend(label(2 "0 -  $1000") label(3 "$1000 - $2000") label(4 "$2000 - $3000") label(5 "$3000 - $4000") label(6 "$4000 - $5000") label(7 "$5000 - $6000")) label(xcoord(longitude) ycoord(latitude) label(F2020))

// Model Selection
test [Wx]lnc = [Wx]lnh = [Wx]lnim = [Wx]lnex = [Wx]lnfd= 0 //SDM-SAR

testnl ([Wx]lnc=-[Spatial]rho*[Main]lnc) ([Wx]lnh=-[Spatial]rho*[Main]lnh) ([Wx]lnim=-[Spatial]rho*[Main]lnim) ([Wx]lnex=-[Spatial]rho*[Main]lnex) ([Wx]lnfd=-[Spatial]rho*[Main]lnfd)

// Most Correct Model is SAR
xsmle lnpy lnc lnh lnex lnim lnfd dbri, wmat(M) model(sar) fe type(ind) effects nolog


//************************************//
// Standard Panel Model Estimation
xtreg dgdp gdp1 expc impc fdic dbri, fe
estimates store fixed
xtreg dgdp gdp1 expc impc fdic dbri, re
estimates store random
hausman fixed random
// 2. Cross Sectional Dependency Test
xtreg lnpy lnc lnh lnex lnim lnfd, fe 
xttest2
xtcsd, pesaran
// Individual
xtcdf lnpy
xtcdf lnc
xtcdf lnh
xtcdf lnex
xtcdf lnim
xtcdf lnfd

// 2. Gen Unit Root Tests
multipurt dgdp gdp1 expc impc fdic, lags(4)

