#tad boundaries

domains <- read.table("./data/GM12878_domain_data_50000.b.txt", header=F)

bounds.GR <- extractBoundaries(domains.mat=domains, 
                               preprocess=FALSE, 
                               CHR=paste0("CHR", c(1:8,10:22)), 
                               resolution=50000)

tadData_oc <- createTADdata(bounds.GR=bounds.GR,
                            resolution=50000,
                            genomicElements.GR=tfbsList,
                            featureType="oc",
                            resampling="none",
                            trainCHR=paste0("CHR",c(1:8,10:21)),
                            predictCHR="CHR22")

all_tad_data <- rbind.data.frame(tadData_oc[[1]], tadData_oc[[2]], stringsAsFactors = FALSE)
all_tad_data <- data.matrix(all_tad_data)
all_tad_data[,1] <- all_tad_data[,1]-1

generator <- function(data, lookback, delay, min_index, max_index,
                      shuffle = FALSE, batch_size = 128, step = 6) {
    if (is.null(max_index))
        max_index <- nrow(data) - delay - 1
    i <- min_index + lookback
    function() {
        if (shuffle) {
            rows <- sample(c((min_index+lookback):max_index), size = batch_size)
        } else {
            if (i + batch_size >= max_index)
                i <<- min_index + lookback
            rows <- c(i:min(i+batch_size-1, max_index))
            i <<- i + length(rows)
        }
        
        samples <- array(0, dim = c(length(rows), 
                                    lookback / step,
                                    dim(data)[[-1]]))
        targets <- array(0, dim = c(length(rows)))
        
        for (j in 1:length(rows)) {
            indices <- seq(rows[[j]] - lookback, rows[[j]] - 1, 
                           length.out = dim(samples)[[2]])
            samples[j,,] <- data[indices,]
            targets[[j]] <- data[rows[[j]] + delay,1]
        }            
        
        list(samples, targets)
    }
}


#setting parameters for the generator function
lookback <- 50 #observations will go back 50 bins (each bin is 50kb in width -> for 50*50000=2500000 bases)
step <- 1 #Observations will be sampled at one data point per bin (50000 bases)
delay <- 10 #Targets will be 10 bins (10*50000=500000 bases) in the future.
batch_size <- 128 #The number of samples per batch
#assigning 3 generator functions for training, validation, and testing
train_gen <- generator(
    all_tad_data,
    lookback = lookback,
    delay = delay,
    min_index = 1,
    max_index = 42561,
    shuffle = FALSE,
    step = step, 
    batch_size = batch_size
)
val_gen = generator(
    all_tad_data,
    lookback = lookback,
    delay = delay,
    min_index = 42562,
    max_index = 52561,
    step = step,
    batch_size = batch_size
)
test_gen <- generator(
    all_tad_data,
    lookback = lookback,
    delay = delay,
    min_index = 52562,
    max_index = NULL,
    step = step,
    batch_size = batch_size
)
# How many steps to draw from val_gen in order to see the entire validation set
val_steps <- (52562 - 42562 - lookback) / batch_size
# How many steps to draw from test_gen in order to see the entire test set
test_steps <- (nrow(all_tad_data) - 52562 - lookback) / batch_size

model <- keras_model_sequential() %>% 
    layer_gru(units = 32, input_shape = list(NULL, dim(all_tad_data)[-1])) %>% 
    layer_dense(units = 1)
model %>% compile(
    optimizer = "rmsprop",
    loss = "binary_crossentropy",
    metrics = c("accuracy")
)
history <- model %>% fit_generator(
    train_gen,
    steps_per_epoch = 500,
    epochs = 10,
    validation_data = val_gen,
    validation_steps = val_steps
)

history_oc_bin <- history

plot(history_oc_bin)

pred_gen_data <- test_gen()
pred_oc <- model %>% predict(pred_gen_data[[1]])
pred_oc <- ifelse(pred_oc < .5,0,1)
table(pred_oc, pred_gen_data[[2]])
