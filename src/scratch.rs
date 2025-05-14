use std::error::Error;
use std::sync::mpsc;
use std::thread;
use std::collections::HashMap;
use std::ops::Deref;
use std::str;
//use std::time::Duration;
use polars::prelude::*;
//use polars_core::prelude::*;
//use polars_io::prelude::*;
//use std::fs::File;
use csv::Writer;
use std::time::Instant;
use fastxgz::fasta_reads;

fn read_into_polars_scratch(file_name: &str) -> PolarsResult<DataFrame> {
    CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file_name.into()))?
        .finish()
}

fn col_mean_scratch(data_in: &DataFrame, col_name: &str) -> f64 {
    let x = data_in.select([col_name]).unwrap().iter().next().unwrap().mean().unwrap();
    x
}

fn get_ref_scratch(file_name: &str) -> Vec<char>{
    let reads = fasta_reads(file_name).expect("The file cannot be opened.");
    let mut x = Vec::new();
    for read in reads {
        x = String::from_utf8(read).unwrap().to_uppercase().chars().collect::<Vec<char>>();
    }
    x
}
fn main() -> Result<(), Box<dyn Error>>{

    let df  = read_into_polars_scratch("somefile.csv")?;
    let x = col_mean_scratch(&df, "methyl_score");
        let num_rows = df.shape().0;
        println!("Input frame has {num_rows} rows and mean score {x}.");
    
        let mut mean_collection = HashMap::new();
        let mut wtr = Writer::from_path("someother_file.csv")?;
        wtr.write_record(&["index","mean"])?;
    
        let chunk_rows = num_rows / 10;
        let mut y = df.slice(0, chunk_rows);
    
        let (tx, rx) = mpsc::channel();
    
        let before = Instant::now();
        for (ind, fr) in y.split_chunks().enumerate() {
            let tx1 = tx.clone();
            thread::spawn(move || {
                tx1.send((ind as u64, col_mean_scratch(&fr, "methyl_score"))).unwrap();
                });
        }
    
        drop(tx);
    
        for (received1, received2) in rx {
            mean_collection.insert(received1, received2);
            wtr.write_record(&[received1.to_string(),received2.to_string()]).unwrap();
            println!("index: {received1}, mean: {received2}");
        }
        println!("Elapsed time: {:.2?}", before.elapsed());
    
    Ok(())
}
