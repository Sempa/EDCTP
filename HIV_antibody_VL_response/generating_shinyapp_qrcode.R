install.packages("rsconnect")
library(rsconnect)

rsconnect::setAccountInfo(name='sempa',
                          token='22479A83F2642CA567F4117568727A0E',
                          secret='LpxqqXJYT1OD+2I4+74Rt88d0lg/Dq5To11nHrHP')
rsconnect::deployApp("HIV_antibody_VL_response")
install.packages("qrcode")
plot(qrcode::qr_code('https://sempa.shinyapps.io/hiv_antibody_vl_response/'))