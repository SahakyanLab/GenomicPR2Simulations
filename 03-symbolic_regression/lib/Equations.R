Equations <- function(data, var){

  # Equations obtained through Eureqa symbolic regression

  # Flag     Format       Description
  # data    <data.frame>  Data set containing rate constants 
  # var     <character>   Mutation rate constant as a variable

  if(!is.character(var)){
    stop("var needs to be a character vector.")
  }
  
  if(var == 'kag'){
    y_pred = 0.348387923987881 + 0.816734963882394*data["ktc"] + 
    0.651414343680873*data["kta"] + 0.609291658916403*data["kgc"] + 
    0.605797402160463*data["kga"] - 0.618906078160872*data["kat"] - 
    0.67329756572978*data["kcg"] - 0.721891543764904*data["kct"]
  } else if(var == 'kat'){
    y_pred = 0.241665477033022 + 0.897500217963755*data["kta"] + 
    0.514297637755254*data["ktg"] + 0.46775568703087*data["ktc"] + 
    0.463842653474543*data["kca"] + 0.393996535254421*data["kga"] - 
    0.479882925548398*data["kag"] - 0.484672348490665*data["kct"] - 
    0.487621790665817*data["kac"] - 0.516076484396038*data["kgt"]
  } else if(var == 'kac'){
    y_pred = 0.22558936599309 + 0.852687563679462*data["kcg"] + 
    0.801456741790903*data["ktg"] + 0.794496499861356*data["kca"] + 
    0.787171624286753*data["kta"] - 0.753938685112867*data["kat"] - 
    0.837621207499993*data["kgc"] - 0.892761983670905*data["kgt"]
  } else if(var == 'kga'){
    y_pred = 0.105020824396609 + 0.736399072631546*data["kag"] + 
    0.736399072631546*data["kat"] + 0.730854831389399*data["kct"] + 
    0.668840893845205*data["kcg"] - 0.643445226498024*data["kgc"] - 
    0.668840893845205*data["ktc"] - 0.738109376474998*data["kta"]
  } else if(var == 'kgt'){
    y_pred = 0.195602186906957 + 0.795944042083753*data["kcg"] + 
    0.762762200816868*data["kta"] + 0.757989214023514*data["ktg"] + 
    0.741388062485789*data["kca"] - 0.709523759605975*data["kat"] - 
    0.758850686673801*data["kgc"] - 0.801757287841475*data["kac"]
  } else if(var == 'kgc'){
    y_pred = 0.0266522601542158 + 0.826978095447326*data["kcg"] + 
    0.489286002530125*data["kag"] + 0.453348892899995*data["kct"] + 
    0.449467306920196*data["kca"] + 0.414417817271229*data["ktg"] - 
    0.385407414591448*data["ktc"] - 0.418343418385765*data["kac"] - 
    0.430961689763262*data["kgt"] - 0.432476148161079*data["kga"]
  } else if(var == 'kta'){
    y_pred = 0.828368985921817*data["kat"] + 0.530245079564161*data["kct"] + 
    0.509302963224534*data["kgt"] + 0.499135884778579*data["kag"] + 
    0.467813865068305*data["kac"] - 0.0460955436591221 - 
    0.42990907018686*data["kca"] - 0.446701251437334*data["ktc"] - 
    0.451430012359511*data["ktg"] - 0.46984114233575*data["kga"]
  } else if(var == 'ktg'){
    y_pred = 0.112574827222794 + 0.88412124636153*data["kgt"] + 
    0.841207273827614*data["kac"] + 0.807204591082419*data["kgc"] + 
    0.794765122731266*data["kat"] - 0.792020815272827*data["kca"] - 
    0.811438663741675*data["kcg"] - 0.835577248635635*data["kta"]
  } else if(var == 'ktc'){
    y_pred = 0.302524932124489 + 0.81778486498004*data["kag"] + 
    0.765641304372685*data["kct"] + 0.620216273339829*data["kat"] + 
    0.605049864248977*data["kcg"] - 0.641586087709269*data["kga"] - 
    0.739907793867365*data["kta"] - 0.757175993792224*data["kgc"]
  } else if(var == 'kca'){
    y_pred = 0.135000854348154 + 0.8899122260307*data["kgt"] + 
    0.826421028836978*data["kgc"] + 0.794875984525169*data["kac"] +
     0.756327474198048*data["kat"] - 0.786605858750487*data["kta"] -
      0.787397685580583*data["ktg"] - 0.840294389850611*data["kcg"]
  } else if(var == 'kcg'){
    y_pred = 0.393304935908064 + 0.764253923040039*data["kgc"] + 
    0.426508013212467*data["kac"] + 0.422715889530294*data["kgt"] + 
    0.415201669034913*data["ktc"] + 0.393304935908064*data["kga"] - 
    0.393304935908064*data["kct"] - 0.417852963826674*data["ktg"] - 
    0.509603310796301*data["kag"] - 0.509603310796301*data["kca"]
  } else if(var == 'kct'){
    y_pred =  0.23960009195787 + 0.787158743182794*data["kgc"] + 
    0.76912313021918*data["kga"] + 0.727489918806558*data["ktc"] + 
    0.695649662821069*data["kta"] - 0.0369776458407348*data["kca"] - 
    0.643942083249785*data["kat"] - 0.71631868285512*data["kcg"] - 
    0.797168947756005*data["kag"]
  }
  
  # true values
  results <- data.frame(
    y_true = data[, var],
    y_pred = y_pred
  )

  return(results)
}