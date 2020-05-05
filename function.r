namerow <- function(x){
  data = x
  rownames(data) = data[,1]
  data = data[,-1]
  if (nrow(data)>4) {
    if(ncol(data)>4){print(data[1:4,1:4])
    }else{
      print(data[1:4,])
    }
  } else {
    if(ncol(data)>4){print(data[,1:4])
    }else{
      print(data)
    }
  }
  return(data)
}

pick <- function(x,y){
  mix=x[rownames(x) %in% rownames(y),]
  if (nrow(mix)>4) {
    if(ncol(mix)>4){print(mix[1:4,1:4])
    }else{
      print(mix[1:4,])
    }
  } else {
    if(ncol(mix)>4){print(mix[,1:4])
    }else{
      print(mix)
    }
  }
  return(mix)
}

readcsv <- function(x){
  data = read.csv(x)
  rownames(data) = data[,1]
  data = data[,-1]
  if (nrow(data)>4) {
    if(ncol(data)>4){print(data[1:4,1:4])
    }else{
      print(data[1:4,])
    }
  } else {
    if(ncol(data)>4){print(data[,1:4])
    }else{
      print(data)
    }
  }
  print(dim(data))
  return(data)
}

see <- function(x){
  if (nrow(x)>4) {
    if(ncol(x)>4){print(x[1:4,1:4])
    }else{
      print(x[1:4,])
    }
  } else {
    if(ncol(x)>4){print(x[,1:4])
    }else{
      print(x)
    }
  }
  print(dim(x))
}