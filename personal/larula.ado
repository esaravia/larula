program larula 
version 13

syntax varlist [if] [in] [, AR(numlist int min=0 max=20) MA(numlist int min=0 max=20) SARIMA(numlist int min=0 max=20)] 


gen n= _n
tsset n

**Graficos
tsline `varlist', name(line, replace) nodraw
ac `varlist', name(ac, replace) nodraw
pac `varlist', name(pac, replace) nodraw
graph combine line ac pac

set more off
tsappend, add(10) 


**1
**NO ESTACIONARIO EN KPSS

kpss `varlist'
scalar r = r(kpss1)
if r >= 0.463 {
	noi di "No estacionario en KPSS"
	tempvar dy
	gen `dy' = D.`varlist'
	dfuller `dy' 
	scalar t = r(p)
		if t <= 0.05 {
					forvalues s=1/`ma'{
					forvalues a=1/`ar'{
					cap arima `varlist', arima (`a',1,`s') 
					estimates store arima
					estat aroots, nograph
					matrix a = r(ar) 
						if a[1,1] <= 0.99 {
						matrix b = r(ma)
							if b[1,1] <= 0.99 {
							**DENTRO DEL CIRCULO UNITARIO
								**HETEROCEDASTICIDAD
								regress `dy'
								estat archlm
								matrix r = r(p)
										if r[1,1] <= 0.05 {
										arch `varlist', arima(`a',1,`s')
										estimates store arch1
										tempvar residual
										predict `residual', res
										corrgram `residual'
										wntestq `residual'
										scalar s = r(p)
										pronosticoc`s', replace
										forecast estimates arch1, names(je1`s'`a', replace)
										forecast solve, begin(_N-10)
										predict mar`s'`a',y dynamic (_N-10) 
										list mar`s'`a' in -11/-1
										ineqdeco mar`s'`a'
										scalar theilar`s'`a' = r(ge1)
											}
											else {
											regress `dy'
											estat archlm
											matrix r = r(p)
											if r[1,1] >= 0.05 {
											tempvar residual
											predict `residual', res
											corrgram `residual'
											wntestq `residual'
											scalar s = r(p)
											if s >= 0.05 {	
											noi di "YAY"
											forecast create pronosticoc`s', replace
											forecast estimates arima, names(je1`s'`a', replace)
											forecast solve, begin(_N-10)
											predict m`s'`a',y dynamic (_N-10) 
											list m`s'`a' in -11/-1
											ineqdeco m`s'`a'
											scalar theil`s'`a' = r(ge1)
											}
						}
					}
				}
			}	
		}
	}
}



if `sarima' >= 2 {	
	tempvar dy
	gen `dy' = D.`varlist'
	dfuller `dy' 
	scalar t = r(p)
		if t <= 0.05 {
					forvalues s=1/`ma'{
					forvalues a=1/`ar' {
					arima `varlist', arima(`a',1,`s') sarima(`a',1,`s',`sarima')
					estimates store arima
					estat aroots, nograph
					matrix a = r(ar) 
						if a[1,1] <= 0.99 {
						matrix b = r(ma)
							if b[1,1] <= 0.99 {
							**DENTRO DEL CIRCULO UNITARIO
								**HETEROCEDASTICIDAD
								regress `dy'
								estat archlm
								matrix r = r(p)
								if r[1,1] >= 0.05 {
									tempvar residual
									predict `residual', res
									**NO TIENE HETEROCEDASTICIDAD
									corrgram `residual'
									wntestq `residual'
									scalar s = r(p)
									if s >= 0.05 {
									noi di "YAY"
									forecast create pronosticoc`s', replace
									forecast estimates arima, names(jes1`s'`a', replace)
									forecast solve, begin(_N-10)
									predict ms`s'`a',y dynamic (_N-10) 
									list ms`s'`a' in -11/-1
									ineqdeco ms`s'`a'
									scalar theils`s'`a' = r(ge1)
								}
						}
					}
				}
			}	
		}
		}
	
}
}
		
**2 
**ESTACIONARIO EN KPSS / NO ESTACIONARIO DICKEY FULLER		
			
kpss `varlist'
scalar r = r(kpss1)
if r <= 0.463 {
**DFULLER
dfuller `varlist'
scalar s = r(p)
	if s >= 0.05 {
		noi di "No estacionario"
		**DIFERENCIAR
		tempvar dy
		gen `dy' = D.`varlist'
		**DFULLER DIFERENCIADO
		dfuller `dy' 
		scalar t = r(p)
		if t <= 0.05 {
			*ESTIMACIONES MA
				forvalues s=1/`ma'{
				forvalues a=1/`ar'{
					cap arima `varlist', arima (`a',1,`s')
					estimates store arima
					estat aroots, nograph
					matrix a = r(ar) 
						if a[1,1] <= 0.99 {
						matrix b = r(ma)
							if b[1,1] <= 0.99 {
							**DENTRO DEL CIRCULO UNITARIO					
							tempvar residual
							predict `residual', res
								**NO AUTOCORRELACION
								regress `residual'
								estat archlm
								matrix r = r(p) 
									if r[1,1] <= 0.05 {
									 arch `varlist', arima(`a',1,`s')
										estimates store arch1
										tempvar residual
										predict `residual', res
										corrgram `residual'
										wntestq `residual'
										scalar s = r(p)
										if s <= 0.05 {
										forecast create pronosticoc`s', replace
										forecast estimates arch1, names(je1`s'`a', replace)
										forecast solve, begin(_N-10)
										predict mar`s'`a',y dynamic (_N-10) 
										list mar`s'`a' in -11/-1
										ineqdeco mar`s'`a'
										scalar theilar`s'`a' = r(ge1)
										}
											else {
											regress `dy'
											estat archlm
											matrix r = r(p)
											if r[1,1] >= 0.05 {
											corrgram `residual'
										    wntestq `residual'
										   scalar s = r(p)
										   if s >= 0.05 {
											noi di "YAY"
											forecast create pronosticoc`s', replace
											forecast estimates arima, names(je2`s'`a', replace)
											forecast solve, begin(_N-10)
											predict m`s'`a', y dynamic (_N-10)
											list m`s'`a' in -11/-1
											ineqdeco m`s'`a'
											scalar theil`s'`a' = r(ge1)
									}
										
									}
								}
							}
				}
			}
		}
	}
}
}



	if `sarima' >= 2 {	
	tempvar dy
	gen `dy' = D.`varlist'
	dfuller `dy' 
	scalar t = r(p)
		if t <= 0.05 {
					forvalues s=1/`ma'{
					forvalues a=1/`ar' {
					arima `varlist', arima(`a',1,`s') sarima(`a',1,`s',`sarima')
					estimates store arima
					estat aroots 
					matrix a = r(ar) 
						if a[1,1] <= 0.99 {
						matrix b = r(ma)
							if b[1,1] <= 0.99 {
							**DENTRO DEL CIRCULO UNITARIO
								**HETEROCEDASTICIDAD
								regress `dy'
								estat archlm
								matrix r = r(p)
								if r[1,1] >= 0.05 {
										tempvar residual
										predict `residual', res
										**NO TIENE HETEROCEDASTICIDAD
										corrgram `residual'
										wntestq `residual'
										scalar s = r(p)
										if s >= 0.05 {
									noi di "YAY"
									forecast create pronosticoc`s', replace
									forecast estimates arima, names(jes1`s'`a', replace)
									forecast solve, begin(_N-10)
									predict ms`s'`a',y dynamic (_N-10) 
									list ms`s'`a' in -11/-1
									ineqdeco ms`s'`a'
									scalar theils`s'`a' = r(ge1)
									}
								}
							}
						}
					}
				}
			}
		}
		}


		
		
**3
**ESTACIONARIO EN KPSS / ESTACIONARIO EN DICKEY FULLER 

		
kpss `varlist'
scalar r = r(kpss1)
if r <= 0.463 {
	dfuller `varlist'
	scalar s = r(p)
	if s <= 0.05 {
	noi di "Estacionario"
		**ESTACIONARIO
			forvalues a=1/`ar' { 
			forvalues s=1/`ma' {
					arima `varlist', arima (`a',0,`s') 
					estimates store arima
					estat aroots 
					matrix a = r(ar) 
						if a[1,1] <= 0.99 {
						matrix b = r(ma)
							if b[1,1] <= 0.99 {
							**DENTRO DEL CIRCULO UNITARIO					
							tempvar residual
							predict `residual', res
								**NO AUTOCORRELACION
								regress `residual'
								estat archlm
								matrix r = r(p) 
									if r[1,1] <= 0.05 {
										arch `varlist', arima(`a',0,`s')
										tempvar residual
										predict `residual', res
										estimates store arch1
										corrgram `residual'
										wntestq `residual'
										scalar s = r(p)
										if s >= 0.05 {
										forecast create pronosticoc`s', replace
										forecast estimates arch1, names(je1`s'`a', replace)
										forecast solve, begin(_N-10)
										predict mar`s'`a',y dynamic (_N-10) 
										list mar`s'`a' in -11/-1
										ineqdeco mar`s'`a'
										scalar theilar`s'`a' = r(ge1)
										}
											else {
											regress `dy'
											estat archlm
											matrix r = r(p)
											if r[1,1] >= 0.05 {
											corrgram `residual'
											wntestq `residual', table
											scalar s = r(p)
											if s >= 0.05 {
											noi di "YAY"
											forecast create pronostico, replace
											forecast estimates arima, names(je3`s'`a', replace)
											forecast solve, begin(_N-10)
											predict m`s'`a', y dynamic (_N-10)
											list m`s'`a' in -11/-1
											ineqdeco m`s'`a'
											scalar theil`s'`a' = r(ge1)
											scalar list theil`s'`a'
									}
										
							}
						}
					}
				}
			}
		}
}
}

	if `sarima' >= 2 {	
	dfuller `varlist'
	scalar t = r(p)
		if t <= 0.05 {
					forvalues s=1/`ma'{
					forvalues a=1/`ar' {
					arima `varlist', arima(`a',0,`s') sarima(`a',0,`s',`sarima')
					estimates store arima
					estat aroots 
					matrix a = r(ar) 
						if a[1,1] <= 0.99 {
						matrix b = r(ma)
							if b[1,1] <= 0.99 {
							**DENTRO DEL CIRCULO UNITARIO
							tempvar residual
							predict `residual', res
							**NO TIENE HETEROCEDASTICIDAD
							corrgram `residual'
							wntestq `residual'
							scalar s = r(p)
							if s >= 0.05 {
								**HETEROCEDASTICIDAD
								regress `dy'
								estat archlm
								matrix r = r(p)
								if r[1,1] >= 0.05 {
									noi di "YAY"
									forecast create pronosticoc`s', replace
									forecast estimates arima, names(jes1`s'`a', replace)
									forecast solve, begin(_N-10)
									predict ms`s'`a',y dynamic (_N-10) 
									list ms`s'`a' in -11/-1
									ineqdeco ms`s'`a'
									scalar theils`s'`a' = r(ge1)
									}
								}
							}
						}
					}
				}
			}
	
}
}


scalar dir


end 


