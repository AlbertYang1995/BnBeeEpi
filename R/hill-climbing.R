
# unified hill climbing implementation (both optimized and by spec).
hill.climbing = function(x, start, whitelist, blacklist, score, extra.args,
    restart, perturb, max.iter, maxp, optimized, debug = FALSE) {

  # cache nodes' labels.
  nodes = names(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # set the iteration counter.
  iter = 1
  # check whether the score is score-equivalent.
  score.equivalence = is.score.equivalent(score, nodes, extra.args)
  # check whether the score is decomposable.
  score.decomposability = is.score.decomposable(score, extra.args)
  # allocate the cache matrix.
  cache.bic = matrix(0, nrow = n.nodes, ncol = n.nodes)
  cache.mit = matrix(0, nrow = n.nodes, ncol = n.nodes)
  # nodes to be updated (all of them in the first iteration).
  updated = seq_len(n.nodes) - 1L
  # set the number of random restarts.
  restart.counter = restart
  # unchanged maximum iter
  unchanged.iter = 4
  # set the population size
  population.size = population.size
  # use fast.iamb result as primary individual
  fast.alpha = fast.alpha
  start = fast.iamb(x, alpha = fast.alpha)

  # convert the blacklist to an adjacency matrix for easy use.
  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)
  else
    blmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)

  # convert the whitelist to an adjacency matrix for easy use.
  if (!is.null(whitelist))
    wlmat = arcs2amat(whitelist, nodes)
  else
    wlmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)

  cat("----------------------------------------------------------------\n")
  cat("* starting from the following network:\n")
  print(start)
  
  # set the metadata of the network; othewise the debugging output is
  # confusing and not nearly as informative.
  start$learning$algo = "hc"
  start$learning$ntests = 0
  start$learning$test = score
  start$learning$args = extra.args
  start$learning$optimized = optimized

  # ====================== begin generating the primary individual ==========================
  # make a copy of fast-iamb result
  start.temp = start
  # the population used to storege all the individuals
  population = array(-1L, dim = c(n.nodes, n.nodes, population.size))  

  # ======================= start generating primary population =========================
  if (narcs(start) == 0) {
    cat("  > The primary network is empty!\n")
    for (i in 1:population.size) {
      cat("  > generate", i, "th individual.\n")
      population[, , i] = arcs2amat(start.temp$arcs, nodes)  
    }
  }
  else {
    arc.size = floor(narcs(start) * 0.3)
    if (arc.size < 1)
      arc.size = 1

    for (i in 1:population.size) {

      random.index = sample(1:narcs(start), size = arc.size)
      fast.amat = amat(start)
      cat("  > generate", i, "th individual.\n")
      for (j in 1:arc.size) {

        k = random.index[j]
        cat("  > delete", arcs(start)[[k, 1]], "-->", arcs(start)[[k, 2]], "\n")
        fast.amat[arcs(start)[[k, 1]], arcs(start)[[k, 2]]] = 0L

      }

      amat(start.temp) = fast.amat
      population[, , i] = arcs2amat(start.temp$arcs, nodes)   

    }
  }
  # ==================== Initializing parameters =========================
  # keep the score of each individual
  population.score = array(-1, population.size)
  # keep the fit score of each individual
  population.fit.score = array(-1, population.size)  
  # used to get the bn. structure 
  network.temp = start
  best.result = start
  # used to record unchanged times
  unchanged.times = array(0, population.size)
  # used to record the individual who has been deleted
  index.record = array(0, population.size)  
  # used to record best individual's score
  best.score = -Inf

  repeat {

    # ============================== picking bee =======================================

    cat("\n  > ============ Begin", iter, "th iteration ============ \n")
    
    # trying to generate second generation

    # change primary individuals randomly, generating second generation    
    for (i in 1:population.size) {
    
      amat(network.temp) = population[, , i]
      nparents = colSums(population[, , i])
      reference.score.bic = per.node.score(network = network.temp, score = "bic",
                                          targets = nodes, extra.args = extra.args, data = x)
      score.bic = sum(reference.score.bic)
      reference.score.mit = per.node.score(network = network.temp, score = "mit",
                                          targets = nodes, extra.args = extra.args, data = x)
      score.mit = sum(reference.score.mit)
      population.score[i] = score.bic + score.mit

      # set up the score cache (BEWARE: in place modification!).
      .Call(call_score_cache_fill,
            nodes = nodes,
            data = x,
            network = network.temp,
            score = "bic",
            extra = extra.args,
            reference = reference.score.bic,
            equivalence = score.equivalence && optimized,
            decomposability = score.decomposability,
            updated = (if (optimized) updated else seq(length(nodes)) - 1L),
            amat = population[, , i],
            cache = cache.bic,
            blmat = blmat,
            debug = debug)

      # set up the score cache (BEWARE: in place modification!).
      .Call(call_score_cache_fill,
            nodes = nodes,
            data = x,
            network = network.temp,
            score = "mit",
            extra = extra.args,
            reference = reference.score.mit,
            equivalence = score.equivalence && optimized,
            decomposability = score.decomposability,
            updated = (if (optimized) updated else seq(length(nodes)) - 1L),
            amat = population[, , i],
            cache = cache.mit,
            blmat = blmat,
            debug = debug)          
      
      # select which arcs should be tested for inclusion in the graph (hybrid
      # learning algorithms should hook the restrict phase here).
      to.be.added = arcs.to.be.added(amat = population[, , i], nodes = nodes,
                                    blacklist = blmat, whitelist = NULL, nparents = nparents,
                                    maxp = maxp, arcs = FALSE)

      # get the best arc addition/removal/reversal.
      bestop = .Call(call_hc_opt_step,
                    amat = population[, , i],
                    nodes = nodes,
                    added = to.be.added,
                    cache_bic = cache.bic,
                    cache_mit = cache.mit,
                    reference_bic = reference.score.bic,
                    reference_mit = reference.score.mit,
                    wlmat = wlmat,
                    blmat = blmat,
                    nparents = nparents,
                    maxp = maxp,
                    debug = debug)
      
      
      if (bestop$op == FALSE) {
        
        cat("----------------------------------------------------------------\n")
        cat("  > for individual", i, ", nothing to do .\n")
        unchanged.times[i] = unchanged.times[i] + 1

        next
        
      }#THEN

      unchanged.times[i] = 0
      
      # update the network structure.
      network.temp = arc.operations(network.temp, from = bestop$from, to = bestop$to,
                                    op = bestop$op, check.cycles = FALSE, check.illegal = FALSE,
                                    update = TRUE, debug = FALSE)
      
      # update the reference score
      population[, , i] = amat(network.temp)
      nparents = colSums(population[, , i])
      reference.score.bic = per.node.score(network = network.temp, score = "bic",
                                          targets = nodes, extra.args = extra.args, data = x)
      score.bic = sum(reference.score.bic)
      reference.score.mit = per.node.score(network = network.temp, score = "mit",
                                          targets = nodes, extra.args = extra.args, data = x)
      score.mit = sum(reference.score.mit)
      population.score[i] = score.bic + score.mit
      
      cat("----------------------------------------------------------------\n")
      if (bestop$op == "set")
        cat("  > add", i, "th individual", bestop$from, "->", bestop$to, ".\n")
      else if (bestop$op == "drop")
        cat("  > drop", i, "th individual", bestop$from, "->", bestop$to, ".\n")
      else
        cat("  > reverse", i, "th individual", bestop$from, "->", bestop$to, ".\n")
      
    }#FOR

    # ========================= caculate fit.score =====================================
    population.fit.score <- sapply(population.score, function(x) {
      result <- 1 / (1 + abs(x))
      return(result)
    })
    score.sum <- sum(population.fit.score)
    population.fit.score <- sapply(population.fit.score, function(x) {
      result <- x / score.sum
      return(result)
    })

    # ================================= observing bee ======================================

    for (chosen.index in 1:population.size) {
      
      if (runif (1, 0, 1) < population.fit.score[chosen.index]) {
        
        cat("----------------------------------------------------------------\n")
        cat("  > trying to update", chosen.index, "th individual.\n")

        amat(network.temp) = population[, , chosen.index]
        nparents = colSums(population[, , chosen.index])
        reference.score.bic = per.node.score(network = network.temp, score = "bic",
                                            targets = nodes, extra.args = extra.args, data = x)
        score.bic = sum(reference.score.bic)
        reference.score.mit = per.node.score(network = network.temp, score = "mit",
                                            targets = nodes, extra.args = extra.args, data = x)
        score.mit = sum(reference.score.mit)
        population.score[chosen.index] = score.bic + score.mit
        
        # set up the score cache (BEWARE: in place modification!).
        .Call(call_score_cache_fill,
              nodes = nodes,
              data = x,
              network = network.temp,
              score = "bic",
              extra = extra.args,
              reference = reference.score.bic,
              equivalence = score.equivalence && optimized,
              decomposability = score.decomposability,
              updated = (if (optimized) updated else seq(length(nodes)) - 1L),
              amat = population[ , , chosen.index],
              cache = cache.bic,
              blmat = blmat,
              debug = debug)

        # set up the score cache (BEWARE: in place modification!).
        .Call(call_score_cache_fill,
              nodes = nodes,
              data = x,
              network = network.temp,
              score = "mit",
              extra = extra.args,
              reference = reference.score.mit,
              equivalence = score.equivalence && optimized,
              decomposability = score.decomposability,
              updated = (if (optimized) updated else seq(length(nodes)) - 1L),
              amat = population[ , , chosen.index],
              cache = cache.mit,
              blmat = blmat,
              debug = debug)

        # select which arcs should be tested for inclusion in the graph (hybrid
        # learning algorithms should hook the restrict phase here).
        to.be.added = arcs.to.be.added(amat = population[, , chosen.index], nodes = nodes,
                                       blacklist = blmat, whitelist = NULL, nparents = nparents,
                                       maxp = maxp, arcs = FALSE)

        # get the best arc addition/removal/reversal.
        bestop = .Call(call_hc_opt_step,
                      amat = population[, , chosen.index],
                      nodes = nodes,
                      added = to.be.added,
                      cache_bic = cache.bic,
                      cache_mit = cache.mit,
                      reference_bic = reference.score.bic,
                      reference_mit = reference.score.mit,
                      wlmat = wlmat,
                      blmat = blmat,
                      nparents = nparents,
                      maxp = maxp,
                      debug = debug)
        

        if (bestop$op == FALSE) {
          
          cat("  > for individual", chosen.index, ", nothing to do .\n")
          unchanged.times[chosen.index] = unchanged.times[chosen.index] + 1
          
          next
          
        }#THEN
        
        unchanged.times[chosen.index] = 0

        # update the network structure.
        network.temp = arc.operations(network.temp, from = bestop$from, to = bestop$to,
                                      op = bestop$op, check.cycles = FALSE, check.illegal = FALSE,
                                      update = TRUE, debug = FALSE)
        
        # update the reference score
        population[, , chosen.index] = amat(network.temp)
        nparents = colSums(population[, , chosen.index])
        reference.score.bic = per.node.score(network = network.temp, score = "bic",
                                            targets = nodes, extra.args = extra.args, data = x)
        score.bic = sum(reference.score.bic)
        reference.score.mit = per.node.score(network = network.temp, score = "mit",
                                            targets = nodes, extra.args = extra.args, data = x)
        score.mit = sum(reference.score.mit)
        population.score[chosen.index] = score.bic + score.mit
        
        if (bestop$op == "set")
          cat("  > add", chosen.index, "th individual", bestop$from, "->", bestop$to, ".\n")
        else if (bestop$op == "drop")
          cat("  > drop", chosen.index, "th individual", bestop$from, "->", bestop$to, ".\n")
        else
          cat("  > reverse", chosen.index, "th individual", bestop$from, "->", bestop$to, ".\n")
        
      }#THEN
      
    }#FOR


    # =============================== detecting bee ======================================
    # record the individual's index who has the maximum unchanged times
    del.index = which.max(unchanged.times)
    
    if (unchanged.times[del.index] > unchanged.iter) {
      
      cat("  > we have to delete", del.index, "th individual and generate a new one.\n")
      cat("  > ============================================= < \n")
      
      # recording which one has been deleted
      index.record[del.index] = 1L

      if (narcs(start) == 0) {
        population[, , del.index] = arcs2amat(start.temp$arcs, nodes)
        cat("  > The primary network is empty!\n")
      }
      else {
        arc.size = floor(narcs(start) * 0.3)
        if (arc.size < 1)
          arc.size = 1
        random.index = sample(1:narcs(start), size = arc.size)
        
        fast.amat = amat(start)

        for (j in 1:arc.size) {

          k = random.index[j]
          cat("  > delete", arcs(start)[[k, 1]], "-->", arcs(start)[[k, 2]], "\n")
          fast.amat[arcs(start)[[k, 1]], arcs(start)[[k, 2]]] = 0L

        }

        amat(start.temp) = fast.amat
        population[, , del.index] = arcs2amat(start.temp$arcs, nodes)
      }
      unchanged.times[del.index] = 0
      amat(network.temp) = population[, , del.index]
      nparents = colSums(population[, , del.index])
      reference.score.bic = per.node.score(network = network.temp, score = "bic",
                                           targets = nodes, extra.args = extra.args, data = x)
      score.bic = sum(reference.score.bic)
      reference.score.mit = per.node.score(network = network.temp, score = "mit",
                                           targets = nodes, extra.args = extra.args, data = x)
      score.mit = sum(reference.score.mit)
      population.score[del.index] = score.bic + score.mit
      
    }#THEN
    
    # record the best individual
    for (i in 1:population.size) {
      
      if (population.score[i] > best.score) {
        
        best.index = i
        best.iter = iter
        best.score = population.score[i]
        # used to record the best solution 
        amat(best.result) = population[, , best.index]
        
      }

    }

    # each individual has been deleted once
    if (sum(index.record) == population.size) break

    # check the current iteration index against the max.iter parameter.
    if (iter >= max.iter) {
      
      if (debug)
        cat("@ stopping at iteration", max.iter, ".\n")
      
      break
      
    }#THEN
    else iter = iter + 1

  }#REPEAT

  cat("  > Finally, the best is", best.index, "th individual in", best.iter,"iteration whose BIC score is: '", best.score, "' .\n")
  cat("  > ============================================= < \n")
  
  print(arcs(best.result))
  return(best.result)

}#HILL.CLIMBING

