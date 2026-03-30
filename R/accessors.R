# ---- R/accessors.R -----------------------------------------------------------
# S4 accessors for GeneRankFit slots
.rfgr_slot_get <- get("slot", envir = asNamespace("methods"))
.rfgr_slot_set <- get("slot<-", envir = asNamespace("methods"))

.rfgr_get <- function(object, name) .rfgr_slot_get(object, name)

.rfgr_set <- function(object, name, value) {
  .rfgr_slot_set(object, name, value = value)
}
.rfgr_params <- function(fit) .rfgr_get(fit, "params")
.rfgr_params_set <- function(fit, value) .rfgr_set(fit, "params", value)

.rfgr_oof <- function(fit) .rfgr_get(fit, "oof")
.rfgr_oof_set <- function(fit, value) .rfgr_set(fit, "oof", value)

.rfgr_imp <- function(fit) .rfgr_get(fit, "imp")
.rfgr_imp_set <- function(fit, value) .rfgr_set(fit, "imp", value)

.rfgr_features <- function(fit) .rfgr_get(fit, "features")
.rfgr_features_set <- function(fit, value) .rfgr_set(fit, "features", value)

.rfgr_final_model <- function(fit) .rfgr_get(fit, "final_model")
.rfgr_final_model_set <- function(fit, value) .rfgr_set(fit, "final_model", value)

.rfgr_calibration <- function(fit) .rfgr_get(fit, "calibration")
.rfgr_calibration_set <- function(fit, value) .rfgr_set(fit, "calibration", value)

.rfgr_model <- function(fit) .rfgr_get(fit, "model")
.rfgr_model_set <- function(fit, value) .rfgr_set(fit, "model", value)