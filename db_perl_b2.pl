    if ($numph > 0) {
					# these are e values used for reagent batch 2 (purified mCP)
                    $cp = 1; # cell path length
                    $ea434=17372;
                    $ea578=94.1;
                    $eb434=2284;
                    $eb578=38676;
                    $Salinity=35;

                    #parse ph

                    $type = hex(substr($phcode,5,2));  # 5-6 char
                    $timesec = hex(substr($phcode,7,8)); #second since 1904/01/01
                      $base_sec = 2082844800;  #second from 1904/01/01 to 1970/01/01
                      $sec_since_19700101 = $timesec - $base_sec;
                      ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($sec_since_19700101);
                      $yr=$year+1900;
                      $mm=$mon+1;
                      $ymdhms = sprintf("%04d/%02d/%02d %02d:%02d:%02d",$yr,$mm,$mday,$hour,$min,$sec);

                    $Temp1 = hex(substr($phcode,15,4));  #temp1
                    $chksum= hex(substr($phcode,463,2)); #the very last 2 characters
                    $Temp2 = hex(substr($phcode,459,4));  #the very last -6 to -3  char.
                    $batt = 0.00366*hex(substr($phcode,455,4));  #the very last -10 to -7  char.
                    $Salinity=hex(substr($phcode,451,4)); # the very last -14 to -11 character -- currently it is 0
                    $Salinity=35;  # don't know why it is set to 35;

                    $Rt1=($Temp1/(4096-$Temp1))*17400;    # Rt1=17400*temp/(4096-temp)
                    $InvT1=0.0010183+0.000241*(log($Rt1))+0.00000015*(log($Rt1))**3;
                    $TempK1=1/$InvT1;
                    $TempC1=$TempK1-273.15;
                    $TempF1=1.8*$TempC1+32;
                    $TempFinal=$TempC1;
                    
                    $Rt2=($Temp2/(4096-$Temp2))*17400;
                    $InvT2=0.0010183+0.000241*(log($Rt2))+0.00000015*(log($Rt2))**3;
                    $TempK2=1/$InvT2;
                    $TempC2=$TempK2-273.15;
                    $TempF2=1.8*$TempC2+32;
                    $T=($TempC1+$TempC2)/2;

                    $TempFinal=($T+($numph-1)*$TempFinal)/$numph;  # calculate mean

					#one value per line, this is the equation used for batch 2 (purified) mCP
                    $pKa=-241.46+0.6367+7085.7/($TempK1)+43.833*(log10($TempK1))-0.080641*($TempK1)
                    -0.3238*($Salinity**0.5)+0.0807*($Salinity)-0.01157*($Salinity**1.5)+0.000694*($Salinity**2);
                    

                    # following value is one per line

                    # Molar absorptivities for reagent batch 2 (purified mCP)
                    $Ea434 = $ea434+(20.612*(24.8-$TempFinal));
                    $Ea578 = $ea578-(1.0177*(24.8-$TempFinal));
                    $Eb434 = $eb434-(6.3863*(24.86-$TempFinal));
                    $Eb578 = $eb578+(66.808*(24.86-$TempFinal));
                    $e1=$Ea578/$Ea434;
                    $e2=$Eb578/$Ea434;
                    $e3=$Eb434/$Ea434;


                    $Ref434A=hex(substr($phcode,19,4));  #0960 = 2400
                    $Sig434A=hex(substr($phcode,23,4));  #042C = 1068
                    $Ref578A=hex(substr($phcode,27,4));  #07A2 = 1954
                    $Sig578A=hex(substr($phcode,31,4));  #0714 = 1814


                    $Ref434B=hex(substr($phcode,35,4));
                    $Sig434B=hex(substr($phcode,39,4));
                    $Ref578B=hex(substr($phcode,43,4));
                    $Sig578B=hex(substr($phcode,47,4));


                    $Ref434C=hex(substr($phcode,51,4));
                    $Sig434C=hex(substr($phcode,55,4));
                    $Ref578C=hex(substr($phcode,59,4));
                    $Sig578C=hex(substr($phcode,63,4));


                    $Ref434D=hex(substr($phcode,67,4));
                    $Sig434D=hex(substr($phcode,71,4));
                    $Ref578D=hex(substr($phcode,75,4));
                    $Sig578D=hex(substr($phcode,79,4));


                    if($Ref434A == 0 || $Ref434B == 0 || $Ref434C == 0 || $Ref434D == 0 ||
                       $Ref578A == 0 || $Ref578B == 0 || $Ref578C== 0 ||  $Ref578D == 0) {
                       $ph_err_flag = -2;

                    }
                    else {
                      $Blank434A=$Sig434A/$Ref434A;
                      $Blank578A=$Sig578A/$Ref578A;
                      $Blank434B=$Sig434B/$Ref434B;
                      $Blank578B=$Sig578B/$Ref578B;
                      $Blank434C=$Sig434C/$Ref434C;
                      $Blank578C=$Sig578C/$Ref578C;
                      $Blank434D=$Sig434D/$Ref434D;
                      $Blank578D=$Sig578D/$Ref578D;

                    }

                    $blank434=($Blank434A+$Blank434B+$Blank434C+$Blank434D)/4;
                    $blank578=($Blank578A+$Blank578B+$Blank578C+$Blank578D)/4;

                    if($blank434 != 0) {
                      $a434blank=-log10($blank434);  #one value per line
                    }
                    else {
                      $a434blank = 'bad';
                    }
                    if($blank578 != 0) {
                      $a578blank=-log10($blank578);
                    }
                    else {
                      $a578blank = 'bad';
                    }
                    $sumy=0;
                    $sumx2=0;
                    $sumy2=0;
                    $sumxy=0;
                    for ($ii=0;$ii<22;$ii++) {
                       $k = 83+16*$ii;
                       $ref434=hex(substr($phcode,$k,4));    #1-4 of 16 chars
                       $i434=hex(substr($phcode,$k+4,4));    #5-8 of 16 chars
                       $ref578=hex(substr($phcode,$k+8,4));  #9-12 of 16 chars
                       $i578=hex(substr($phcode,$k+12,4));    #12-16 of 16 chars

                       if ($i434 == 0 || $i578 == 0 || $a578blank == 'bad' || $a434blank == 'bad') {
                         $ph_err_flag = -2;
                       }
                       else {
                         if($ref434 == 0 || $ref578 == 0) {
                            $ph_err_flag = -2;
                         }
                         else {
                           # Absorbance
                           $a434 = -log10($i434/$ref434);   #22 values
                           $a578 = -log10($i578/$ref578);   #22 values
                           $abs434=$a434-$a434blank;
                           $abs578=$a578-$a578blank;
                           if ($abs434 != 0) {
                             $R = $abs578 / $abs434;
                             $V1 = $R - $e1;
                             $V2 = $e2 - $R * $e3;
                           }
                           else {
                             $ph_err_flag = -2;
                           }
                        }
							#I changed these equations; you need to use the absorbances, not the ratios
                         $HI=(($abs434*$Eb578)-($abs578*$Eb434))/(($Ea434*$Eb578)-($Eb434*$Ea578)); #22 values
                         $I=(($abs578*$Ea434)-($abs434*$Ea578))/(($Ea434*$Eb578)-($Eb434*$Ea578));  #22 values

                         if ($ii == 13) {
                            $xbar=0;
                            $ybar=0;
                         }
                         # I'm not sure exactly what this does, but I think it does a regression of points 15
                         # to 22; what it SHOULD do instead is go through points 5 through 22, and find the region
                         # with 5 continuous points which produce the best R^2, and have 0.00004<IndConc<0.00008.
                         
                         if ($ii>13 && $ii < 22) { #cutoff=15 to 22
                            # Use data points that are in linear region
                            $IndConc=$HI+$I;  # only use 15:22; but will get all there
                            $pointpH = $pKa+log10($V1/$V2);
                            # Use data points that are in linear region
                            $sumx += $IndConc;
                            $sumy += $pointpH;
                            $sumxy += $IndConc * $pointpH;
                            $sumx2 += $IndConc **2;
                            $sumy2 += $pointpH **2;
                            $xbar = ($IndConc+($ii-14)*$xbar)/($ii-13);
                            $ybar = ($pointpH+($ii-14)*$ybar)/($ii-13);
                         }
                       }
                    }
                    if ($ph_err_flag != -2) {
                      $s2=8;
                      $sumxx2=$sumx*$sumx;   # single value
                      $sumyy2=$sumy*$sumy;

                      $ssxy=$sumxy-($sumx*$sumy)/$s2;
                      $ssx=$sumx2-($sumxx2/$s2);
                      $ssy=$sumy2-($sumyy2/$s2);

                      $slope=$ssxy/$ssx;
                      $finalpH=$ybar-($slope*$xbar);
                      $r2=(($ssxy**2)/($ssx*$ssy));
                      $temp = sprintf("%7.4f",$TempFinal);
                      $ph = sprintf("%7.4f",$finalpH);

                      # Reality check computed pH
                      if ($finalpH >= 7.0 && $finalpH <= 8.9) {  # or 7.5 to 8.4
                        $ph_err_flag = 0;
                        $thephflag=0;
                      }
                      else {
                        $ph_err_flag = 11; # bad data
                        $thephflag = 11;

                      }
                      if ($ph_flag == 12) { #if sensor has problem, then we set flag=11;
                        $ph_err_flag = 12;
                        $thephflag=11;
                      }
                    }
                    $insert_ph = "insert into ph_$tblid (fileid,hdrtime,phtemp,ph,batt,phtime,flag,phstr,phflag) values "."($curr_ID,'$phhdrtime',$temp,$ph,$batt,'$ymdhms',$ph_flag,'$phcode',$thephflag)";
                  }
              }