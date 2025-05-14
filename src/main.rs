mod scratch;

use std::error::Error;
use std::sync::mpsc;
use std::thread;
use std::collections::HashMap;
use std::ops::Deref;
use std::str;
//use std::time::Duration;
use polars::prelude::*;
//use polars_core::prelude::*;
use polars_io::prelude::*;
use std::fs::File;
use csv::Writer;
use std::time::Instant;
use fastxgz::fasta_reads;


//#[cfg(feature = "parquet")]

const INPUT_FILE: &str = "/Users/sveta/xella_consulting/data/chr6_ma_calls_p012.csv";
const PARQUET_FILE: &str = "TXHB102-0016.all.parquet";
const DATA_PATH: &str = "/Users/sveta/xella_consulting/data/rust_project/";
const OUTPUT_FILE: &str = "output.csv";
const CHR: &str = "chr6";
const REF_FILE: &str = "chr6.fa.gz";
const COLS: [&str; 8]  = ["#ct", "st", "en", "strand", "fiber", "fiber_length", "ref_m6a", "m6a_qual"];
const WIN_OFFSET: i64 = 50000;
const INTERVAL: i64 = 1000;
const MIN_DEPTH: usize = 8;
const MIN_SCORE: f64 = 120.;


fn read_into_polars(file_name: &str) -> PolarsResult<DataFrame> {
    CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file_name.into()))?
        .finish()
}

fn read_into_polars2(file_name: &str) -> PolarsResult<DataFrame> {
    LazyFrame::scan_parquet(file_name, ScanArgsParquet::default())?
        .select([cols(COLS)])
        .filter(col("#ct").eq(lit(CHR)))
        .collect()
    //let mut file = File::open([DATA_PATH,PARQUET_FILE].join("")).unwrap();
    //ParquetReader::new(&mut file).finish()
}

fn col_mean(data_in: &DataFrame, col_name: &str) -> f64 {
    let x = data_in.select([col_name]).unwrap().iter().next().unwrap().mean().unwrap();
    x
}

fn get_ref(file_name: &str) -> Vec<char>{
    let reads = fasta_reads(file_name).expect("The file cannot be opened.");
    let mut x = Vec::new();
    for read in reads {
        x = String::from_utf8(read).unwrap().to_uppercase().chars().collect::<Vec<char>>();
    }
    x
}

fn process_region(region_data: &DataFrame, region_start: i64, region_end: i64) -> PolarsResult<DataFrame> {
    let region_depth = region_data.shape().0;
    let schema = Schema::from_iter(vec![
        Field::new("coords".into(), DataType::Int64),
        Field::new("norm_score".into(), DataType::Float64),
        Field::new("strand".into(), DataType::String),
    ]);
    let mut result_df = DataFrame::empty_with_schema(&schema);

    if region_depth < MIN_DEPTH {
        return Ok(result_df);
    }
    let num_steps = (region_end - region_start) / INTERVAL;
    let mut window_end = region_start;
    for i in 0..num_steps {
        let window_start = window_end;
        window_end = window_start + INTERVAL;
        let df_loc = region_data.clone().lazy()
            .filter(col("st").gt(lit(window_start - WIN_OFFSET))
                .and(col("en").lt(lit(window_end + WIN_OFFSET))))
            .filter(col("ref_m6a").eq(lit(".")).not()
                .and(col("m6a_qual").eq(lit(".")).not()))
            .collect().unwrap();
        let window_depth = df_loc.shape().0;
        if window_depth <= MIN_DEPTH {
            continue;
        }
        let df_loc_plus = df_loc
            .clone().lazy().filter(col("strand").eq(lit("+"))).collect()?;
        let df_loc_minus = df_loc
            .clone().lazy().filter(col("strand").eq(lit("-"))).collect()?;
        let window_depth_plus = df_loc_plus.shape().0;
        if window_depth_plus > MIN_DEPTH {
            let cr = df_loc_plus.clone().lazy()
                .select([col("ref_m6a")]).collect()?;
            let qr = df_loc_plus.clone().lazy()
                .select([col("m6a_qual")]).collect()?;
            let mut y: Vec<i64> = Vec::new();
            let mut z: Vec<i64> = Vec::new();
            let mut norm_depth = 0.01;
            for i in 0..cr.height(){
                let yy = cr.get_row(i).iter().next()
                    .unwrap().clone().0.iter()
                    .map(|x| x.to_string()).collect::<Vec<String>>()
                    .pop().unwrap().replace("\"","").strip_suffix(",").unwrap().to_string();
                let yyy: Vec<i64> = yy.split(",").map(|x| x.parse::<i64>().unwrap()).collect();
                norm_depth = if *yyy.iter().min().unwrap() > window_end || 
                    *yyy.iter().max().unwrap() < window_start { norm_depth
                } else {norm_depth + 1.};
                let zz = qr.get_row(i).iter().next()
                    .unwrap().clone().0.iter()
                    .map(|x| x.to_string()).collect::<Vec<String>>()
                    .pop().unwrap().replace("\"","").strip_suffix(",").unwrap().to_string();
                let zzz: Vec<i64> = zz.split(",").map(|x| x.parse::<i64>().unwrap()).collect();
                y.extend(&yyy);
                z.extend(&zzz);
            }
            let norm_factor = 1.0 / (norm_depth);
            let coord_col = Column::new("coords".into(), y);
            let qual_col = Column::new("score".into(), z);
            let m6a_df_full = DataFrame::new(vec![coord_col, qual_col])?;
            let m6a_df = m6a_df_full.clone().lazy().group_by([col("coords")])
                .agg([(col("score")*lit(norm_factor)).sum().alias("norm_score")])
                .filter(col("norm_score").gt(lit(MIN_SCORE)))
                .filter(col("coords").gt(lit(window_start))
                            .and(col("coords").lt(lit(window_end))))
                .with_column(lit("+").alias("strand"))
                .collect()?;
            result_df = concat([result_df.clone().lazy(), m6a_df.clone().lazy()],
                UnionArgs::default(),)?.collect()?;
        }
        let window_depth_minus = df_loc_minus.shape().0;
        if window_depth_minus > MIN_DEPTH {
            let cr = df_loc_minus.clone().lazy()
                .select([col("ref_m6a")]).collect()?;
            let qr = df_loc_minus.clone().lazy()
                .select([col("m6a_qual")]).collect()?;
            let mut y: Vec<i64> = Vec::new();
            let mut z: Vec<i64> = Vec::new();
            let mut norm_depth = 0.01;
            for j in 0..cr.height(){
                let yy = cr.get_row(j).iter().next()
                    .unwrap().clone().0.iter()
                    .map(|x| x.to_string()).collect::<Vec<String>>()
                    .pop().unwrap().replace("\"","").strip_suffix(",").unwrap().to_string();
                let yyy: Vec<i64> = yy.split(",").map(|x| x.parse::<i64>().unwrap()).collect();
                norm_depth = if *yyy.iter().min().unwrap() > window_end ||
                    *yyy.iter().max().unwrap() < window_start { norm_depth
                } else {norm_depth + 1.};
                let zz = qr.get_row(j).iter().next()
                    .unwrap().clone().0.iter()
                    .map(|x| x.to_string()).collect::<Vec<String>>()
                    .pop().unwrap().replace("\"","").strip_suffix(",").unwrap().to_string();
                let zzz: Vec<i64> = zz.split(",").map(|x| x.parse::<i64>().unwrap()).collect();
                y.extend(&yyy);
                z.extend(&zzz);
            }
            let norm_factor = 1.0 / (norm_depth);
            let coord_col = Column::new("coords".into(), y);
            let qual_col = Column::new("score".into(), z);
            let m6a_df_full = DataFrame::new(vec![coord_col, qual_col])?;
            let m6a_df = m6a_df_full.clone().lazy().group_by([col("coords")])
                .agg([(col("score")*lit(norm_factor)).sum().alias("norm_score")])
                .filter(col("norm_score").gt(lit(MIN_SCORE)))
                .filter(col("coords").gt(lit(window_start))
                    .and(col("coords").lt(lit(window_end))))
                .with_column(lit("-").alias("strand"))
                .collect()?;
            result_df = concat([result_df.clone().lazy(), m6a_df.clone().lazy()],
                               UnionArgs::default(),)?.collect()?;
        }
    }
    Ok(result_df)
}

fn main() -> Result<(), Box<dyn Error>>{

    let df  = read_into_polars2(&[DATA_PATH,PARQUET_FILE].concat())?;
    let num_intervals: i64 = 100;
    let mut r_end: i64 = 13075000;
    let mut r_begin: i64;
    let num_threads: usize = 2;
    let mut file = File::create(&[DATA_PATH,OUTPUT_FILE].concat()).expect("could not create file");

    let schema = Schema::from_iter(vec![
        Field::new("coords".into(), DataType::Int64),
        Field::new("norm_score".into(), DataType::Float64),
        Field::new("strand".into(), DataType::String),
    ]);
    let mut final_result_df = DataFrame::empty_with_schema(&schema);
    
    let (tx, rx) = mpsc::channel();

    let before = Instant::now();

    for ind in 0..num_threads {
        r_begin = r_end;
        r_end = r_begin + num_intervals * INTERVAL;

        let df_loc = df.clone().lazy()
            .filter(col("st").gt(lit(r_begin - WIN_OFFSET))
                .and(col("en").lt(lit(r_end + WIN_OFFSET))))
            .collect().unwrap();
        let tx1 = tx.clone();
        thread::spawn(move || {
        tx1.send((ind as u64, process_region(&df_loc, r_begin, r_end))).unwrap();
        // let r_df = process_region(&df_loc, r_begin, r_end)?;
        });
    }
    
    drop(tx);

    for (received1, received2) in rx {
        println!("Received result from thread {received1}");
        let result_df = received2?.clone();
        if result_df.shape().0 > 0 {
            final_result_df = concat([final_result_df.clone().lazy(), result_df.clone().lazy()],
                                     UnionArgs::default(), )?.collect()?;
        }
    }
    println!("Elapsed time: {:.2?}", before.elapsed());
    println!("Result: {:?}", final_result_df.head(None));
    
    CsvWriter::new(&mut file).include_header(true).with_separator(b',')
        .finish(&mut final_result_df);

    Ok(())
}
