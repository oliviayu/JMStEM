#' Update columns grouped by id
#'
update_columns <- function(new_dat, id, ..., list_input=NULL){
  data_list <- c(list(...), list_input)
  cols <- setdiff(names(new_dat), id)
  res <- lapply(data_list, function(dati){
    dplyr::select(dati, -tidyselect::any_of(cols)) %>%
      dplyr::left_join(new_dat, by = id)
  })

  res
}


calculate_fixeff <- function(dat1, id=NULL, equations, pars, names){
  uniqueID <- unique(dat1[[id]])

  eta_str <- group_by_at(dat1, id) %>%
    group_map( function(dt, ...){
      x <- lapply(equations, function(md){
        model.matrix(md, data= dt[1,])
      }) %>% Matrix::bdiag() %>% as.matrix()
      res <- x %*% matrix(pars)
      as.numeric(res)
    }) %>% do.call(rbind, .) %>% as.data.frame()
  names(eta_str) <- names
  eta_str[, id] <- uniqueID
  eta_str
}
