use std::f32::INFINITY;
use bitlab::*;

use std::fs;
use fs::File;
use std::io::Read;

pub fn read_to_vec(filename: &String) -> Vec<u8> {
  let mut f = File::open(&filename).expect("File doesn't exist.");
  let metadata = fs::metadata(&filename).expect("Unable to read metadata.");
  let mut buffer = vec![0; metadata.len() as usize];
  f.read(&mut buffer).expect("Buffer overflow.");
  buffer
}


pub struct DimRanges([f32; 2], [f32; 2], [f32; 2]);


pub struct Transform {}
impl Transform {
  pub fn to_cartesian([r, p, t]: [f32; 3]) -> [f32; 3] {
    [
      r*p.sin()*t.sin(),
      r*p.cos(),
      r*p.sin()*t.cos(),
    ]
  }
  pub fn to_noised_spherical() {

  }
  pub fn to_noised_cartesian() {
    println!("to_noised_cartesian");
  }

  pub fn from_custom_triplet() {
    println!("from_custom_triplet");
  }
}

pub struct Header {}
impl Header {
  pub fn read_u3(vec: &Vec<u8>) -> [u8; 5] {
    // We know there are 5 in the header spec, so hard-coding:
    let mut u3: [u8; 5] = [0; 5];
    for i in 0..5 {
      u3[i as usize] = vec.get_u8(0, i*3, 3).unwrap();
    }
    u3
  }
}


pub fn fast_map_triplets(mut coords: &Vec<f32>, mapper: &dyn Fn([f32; 3]) -> [f32; 3]) -> Vec<f32> {
  vec![1.,2.]
}

pub fn get_dimension_ranges(nums: &Vec<f32>) -> [f32; 3] {
  let mut max = [-INFINITY, -INFINITY, -INFINITY];
  let mut min = [INFINITY, INFINITY, INFINITY];
  min
}


#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn to_cartesian_works() {
    let [r, p, t] = [1., 2., 3.];
    assert_eq!(
      Transform::to_cartesian([r,p,t]),
      [
        r*p.sin()*t.sin(),
        r*p.cos(),
        r*p.sin()*t.cos(),
      ],
    );
  }
  #[test]
  fn dimension_ranges() {
    assert_eq!(
      get_dimension_ranges(&vec![1.,2.,3.]),
      [INFINITY, INFINITY, INFINITY],
    )
  }
  #[test]
  fn bitlab_works() {
    // use bitlab::*;
    let a: i8 = -33; // 0b1101_1111;
    let b = a.get_u8(1,3).unwrap();
    assert_eq!(b, 5);
  }
  #[test]
  fn read_to_vec_works() {
    read_to_vec(&String::from("big.bin"));
  }

  #[test]
  fn read_3bit_sanity_check() {
    let mut u3: [u8; 5] = [0; 5];
    for i in 0..5 {
      let vec = read_to_vec(&String::from("big.bin"));
      let num = vec.get_u8(0, i*3, 3).unwrap();
      u3[i as usize] = num;
      // let vec = read_to_vec(&String::from("big.bin"));
      // let mb1 = vec.get_u8(0, i, 3).unwrap();
      println!("{} #{}", num, i*3);
      println!(" ");
    }
    // for i in 0..3 {
    //   let vec = read_to_vec(&String::from("tiny.bin"));
    //   let mb1_2 = vec.get_u8(1, i*3, 3).unwrap();

    //   // let vec = read_to_vec(&String::from("big.bin"));
    //   // let mb1 = vec.get_u8(0, i, 3).unwrap();
    //   println!("{} #{}", mb1_2+1, i*3);
    //   println!(" ");
    // }


  }
}