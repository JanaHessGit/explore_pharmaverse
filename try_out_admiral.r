# admiral - ADaM dataset construction
# metacore - Stores and validates dataset metadata (variable definitions, codelists) from specs
# xortr - Applies SAS XPT metadata (types, labels, lengths, formats) and exports submission-ready .xpt transport files
# rtables - Framework for building complex nested tables with flexible layouts


library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(stringr)
library(tibble)
library(pharmaversesdtm)
library(admiral)
library(languageserver)
library(admiralonco)
library(metacore)
library(metatools)
library(xportr)
library(haven)


# Read in SDTM datasets
dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
vs <- pharmaversesdtm::vs
admiral_adsl <- admiral::admiral_adsl


ex_ext <- ex %>% derive_vars_dtm(dtc = EXSTDTC, new_vars_prefix = "EXST")
colnames(ex)
colnames(ex_ext)
# columns EXSTDTM and EXSTTMF added: imputed datetime *DTM and the date/time imputation flag 

vs <- vs %>% filter(USUBJID %in% c("01-701-1015", "01-701-1023", "01-703-1086", "01-703-1096", "01-707-1037", "01-716-1024") &
                    VSTESTCD %in% c("SYSBP", "DIABP") &
                    VSPOS == "SUPINE"
                    )

adsl <- admiral_adsl %>% select(-TRTSDTM, -TRTSTMF)

advs <- vs %>% mutate(PARAM = VSTEST, PARAMCD = VSTESTCD, AVAL = VSSTRESN, AVALU = VSORRESU, AVISIT = VISIT, AVISITN = VISITNUM)



# add the treatment start datetime and corresponding imputation flag (EXSTTMF) to adsl by identifying the first record in ex with a non-missing Exposure Start Datetime (EXSTDTM) when sorting by EXSTDM and EXSEQ.
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  )


# derive the Mean Arterial Pressure (MAP) from the Systolic and Diastolic values of the blood pressure
advs <- advs %>%
    derive_param_computed(
        by_vars = exprs(USUBJID, AVISIT, AVISITN), # groups -> one calculation per subj and visit
        parameters = c("SYSBP", "DIABP"),
        set_values_to = exprs(
            AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
            PARAMCD = "MAP",
            PARAM = "Mean Arterial Pressure (mmHg)",
            AVALU = "mmHg"
        )
    )


# Add the date of the final lab visit to ADSL
adsl <- dm %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(FINLABDT = convert_dtc_to_dt(DSSTDTC)),
    filter_add = DSDECOD == "FINAL LAB VISIT"
  )
# Computations can also be used inside a mutate() statement (or similar), so they can be leveraged when your variable derivations do not require calls to specific admiral functions:
  adsl <- adsl %>% mutate(RFSTDT = convert_dtc_to_dt(RFSTDTC))




advs_lastmap <- advs %>% 
  filter(PARAMCD == "MAP") %>%
  #  extract the most recent MAP records
  filter_extreme(
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN, PARAMCD),
    mode = "last"
  )



# expressions
get_admiral_option("subject_keys")
adcm1 <- data.frame(STUDYID = "XXX", USUBJID = "XXX-1", CMTRT = "ASPIRIN")
adcm1

# This doesn't work as we are not unquoting the subject keys
# adcm1 %>% select(get_admiral_option("subject_keys"))
# This works because we are unquoting the subject keys
adcm1 %>% select(!!!get_admiral_option("subject_keys"))



vs1 <- tribble(
  ~USUBJID, ~VSTESTCD, ~VISIT, ~VSSTRESN, ~VSSTRESU, ~VSDTC,
  "01-1301", "WEIGHT", "SCREENING", 82.1, "kg", "2013-08-29",
  "01-1301", "WEIGHT", "WEEK 2", 81.19, "kg", "2013-09-15",
  "01-1301", "WEIGHT", "WEEK 4", 82.56, "kg", "2013-09-24",
  "01-1302", "BMI", "SCREENING", 20.1, "kg/m2", "2013-08-29",
  "01-1302", "BMI", "WEEK 2", 20.2, "kg/m2", "2013-09-15",
  "01-1302", "BMI", "WEEK 4", 19.9, "kg/m2", "2013-09-24"
)

dm1 <- tribble(
  ~USUBJID, ~AGE,
  "01-1301", 18
)


my_expression <- expr(VSTESTCD == "WEIGHT" & VISIT == "SCREENING")

derive_vars_merged(
  dm1,
  dataset_add = select(vs1, USUBJID, VSTESTCD, VISIT),
  by_vars = exprs(USUBJID),
  filter_add = !!my_expression
)



# Derive birth date from BRTHDTC
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "BRTH",
    dtc = BRTHDTC
  )

#decide which reference date to use
colnames(adsl)[grepl("DTC", colnames(adsl))]

any(is.na(adsl$RFSTDTC))
any(is.na(adsl$RFXENDTCDT))
any(is.na(adsl$RFICDTC))
any(is.na(adsl$DMDTC))




 adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "DMDTC",
    dtc = DMDTC
  )
  
colnames(adsl)


adsl <- adsl %>%
  derive_vars_aage(
    start_date = BRTHDT,
    end_date = DMDTCDT
  )


#compare if that was the ref date used
any(is.na(adsl$AAGE))

cmp <- adsl %>%
  transmute(
    USUBJID,
    AGE,
    AAGE,
    age_diff = AGE - AAGE,
    matches = AGE == AAGE
  )
  
all(cmp$matches)
  
  

# understand idea of metacore

load(metacore_example("pilot_ADaM.rda"))
metacore
metacore <- metacore %>% select_dataset("ADSL")
ds_list <- list(DM = read_xpt(metatools_example("dm.xpt")))
ds_list
lapply(ds_list$DM, attr, which = "label")

# Pull in columns from DM to be in ADSL
adsl <- build_from_derived(metacore, ds_list, predecessor_only = FALSE)  

# Now we can make some new columns and converting columns into factors
adsl %>%
   select(USUBJID, AGE, ETHNIC) %>% 
   create_cat_var(metacore, AGE, AGEGR1, AGEGR1N) %>% #Add an AGEGR1 and AGEGR1N column
   convert_var_to_fct(metacore, ETHNIC) # Change ETHNIC to as factor 
