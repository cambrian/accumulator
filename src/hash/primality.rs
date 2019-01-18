use num::bigint::BigInt;
use num::bigint::BigUint;
use num_traits::cast::ToPrimitive;

#[allow(dead_code)]
const SMALL_PRIMES: [u64; 50] = [
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
  101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
  197, 199, 211, 223, 227, 229,
];

// Baillie-PSW probabilistic primality test:
// 1. Filter composites with small divisors.
// 2. Do Miller-Rabin base 2.
// 3. Filter squares.
// 4. Do Lucas.
#[allow(dead_code)]
pub fn is_prob_prime(n: &BigUint) -> bool {
  !has_small_prime_factor(n) && passes_miller_rabin_base_2(n) && !is_square(n) && passes_lucas(n)
}

#[allow(dead_code)]
fn has_small_prime_factor(n: &BigUint) -> bool {
  for &divisor in SMALL_PRIMES.iter() {
    let divisor = &BigUint::from(divisor);
    if divisor > n {
      break;
    }
    if n % divisor == BigUint::from(0u64) {
      return true;
    }
  }
  false
}

#[allow(dead_code)]
fn passes_miller_rabin_base_2(n: &BigUint) -> bool {
  let one = BigUint::from(1u64);
  let two = BigUint::from(2u64);

  // write n-1 = 2^r * d
  let mut d = n - 1u64;
  let mut r = 0;
  while &d % 2u64 == BigUint::from(0u64) {
    d /= 2u64;
    r += 1;
  }
  // println!("{} = 2^{} * {}", n, r, d);
  let mut x = two.modpow(&d, n);
  if x == one || x == n - &one {
    return true;
  }
  for _ in 0..(r - 1) {
    x = x.modpow(&two, n);
    if x == one {
      return false;
    }
    if x == n - &one {
      return true;
    }
  }
  false
}

fn is_square(n: &BigUint) -> bool {
  // Step 1
  let zero = BigUint::from(0 as u8);
  let one = BigUint::from(1 as u8);
  if n & BigUint::from(2 as u8) != zero
    || n & BigUint::from(7 as u8) == BigUint::from(5 as u8)
    || n & BigUint::from(11 as u8) == BigUint::from(8 as u8)
  {
    return false;
  }
  // Maybe unneccessary
  if *n == zero {
    return true;
  }

  let bad255 = [
    false, false, true, true, false, true, true, true, true, false, true, true, true, true, true,
    false, false, true, true, false, true, false, true, true, true, false, true, true, true, true,
    false, true, true, true, false, true, false, true, true, true, true, true, true, true, true,
    true, true, true, true, false, true, false, true, true, true, false, true, true, true, true,
    false, true, true, true, false, true, false, true, true, false, false, true, true, true, true,
    true, false, true, true, true, true, false, true, true, false, false, true, true, true, true,
    true, true, true, true, false, true, true, true, true, true, false, true, true, true, true,
    true, false, true, true, true, true, false, true, true, true, false, true, true, true, true,
    false, false, true, true, true, true, true, true, true, true, true, true, true, true, true,
    false, false, true, true, true, true, true, true, true, false, false, true, true, true, true,
    true, false, true, true, false, true, true, true, true, true, true, true, true, true, true,
    true, false, true, true, false, true, false, true, true, false, true, true, true, true, true,
    true, true, true, true, true, true, false, true, true, false, true, true, true, true, true,
    false, false, true, true, true, true, true, true, true, false, false, true, true, true, true,
    true, true, true, true, true, true, true, true, true, false, false, true, true, true, true,
    false, true, true, true, false, true, true, true, true, false, true, true, true, true, true,
    false, true, true, true, true, true, false, true, true, true, true, true, true, true, true,
    false, false, true, true, false, true, true, true, true, false, true, true, true, true, true,
    false, false, true, true, false, true, false, true, true, true, false, true, true, true, true,
    false, true, true, true, false, true, false, true, true, true, true, true, true, true, true,
    true, true, true, true, false, true, false, true, true, true, false, true, true, true, true,
    false, true, true, true, false, true, false, true, true, false, false, true, true, true, true,
    true, false, true, true, true, true, false, true, true, false, false, true, true, true, true,
    true, true, true, true, false, true, true, true, true, true, false, true, true, true, true,
    true, false, true, true, true, true, false, true, true, true, false, true, true, true, true,
    false, false, true, true, true, true, true, true, true, true, true, true, true, true, true,
    false, false, true, true, true, true, true, true, true, false, false, true, true, true, true,
    true, false, true, true, false, true, true, true, true, true, true, true, true, true, true,
    true, false, true, true, false, true, false, true, true, false, true, true, true, true, true,
    true, true, true, true, true, true, false, true, true, false, true, true, true, true, true,
    false, false, true, true, true, true, true, true, true, false, false, true, true, true, true,
    true, true, true, true, true, true, true, true, true, false, false, true, true, true, true,
    false, true, true, true, false, true, true, true, true, false, true, true, true, true, true,
    false, true, true, true, true, true, false, true, true, true, true, true, true, true, true,
    false, false,
  ];

  let start = [
    1, 3, 1769, 5, 1937, 1741, 7, 1451, 479, 157, 9, 91, 945, 659, 1817, 11, 1983, 707, 1321, 1211,
    1071, 13, 1479, 405, 415, 1501, 1609, 741, 15, 339, 1703, 203, 129, 1411, 873, 1669, 17, 1715,
    1145, 1835, 351, 1251, 887, 1573, 975, 19, 1127, 395, 1855, 1981, 425, 453, 1105, 653, 327, 21,
    287, 93, 713, 1691, 1935, 301, 551, 587, 257, 1277, 23, 763, 1903, 1075, 1799, 1877, 223, 1437,
    1783, 859, 1201, 621, 25, 779, 1727, 573, 471, 1979, 815, 1293, 825, 363, 159, 1315, 183, 27,
    241, 941, 601, 971, 385, 131, 919, 901, 273, 435, 647, 1493, 95, 29, 1417, 805, 719, 1261,
    1177, 1163, 1599, 835, 1367, 315, 1361, 1933, 1977, 747, 31, 1373, 1079, 1637, 1679, 1581,
    1753, 1355, 513, 1539, 1815, 1531, 1647, 205, 505, 1109, 33, 1379, 521, 1627, 1457, 1901, 1767,
    1547, 1471, 1853, 1833, 1349, 559, 1523, 967, 1131, 97, 35, 1975, 795, 497, 1875, 1191, 1739,
    641, 1149, 1385, 133, 529, 845, 1657, 725, 161, 1309, 375, 37, 463, 1555, 615, 1931, 1343, 445,
    937, 1083, 1617, 883, 185, 1515, 225, 1443, 1225, 869, 1423, 1235, 39, 1973, 769, 259, 489,
    1797, 1391, 1485, 1287, 341, 289, 99, 1271, 1701, 1713, 915, 537, 1781, 1215, 963, 41, 581,
    303, 243, 1337, 1899, 353, 1245, 329, 1563, 753, 595, 1113, 1589, 897, 1667, 407, 635, 785,
    1971, 135, 43, 417, 1507, 1929, 731, 207, 275, 1689, 1397, 1087, 1725, 855, 1851, 1873, 397,
    1607, 1813, 481, 163, 567, 101, 1167, 45, 1831, 1205, 1025, 1021, 1303, 1029, 1135, 1331, 1017,
    427, 545, 1181, 1033, 933, 1969, 365, 1255, 1013, 959, 317, 1751, 187, 47, 1037, 455, 1429,
    609, 1571, 1463, 1765, 1009, 685, 679, 821, 1153, 387, 1897, 1403, 1041, 691, 1927, 811, 673,
    227, 137, 1499, 49, 1005, 103, 629, 831, 1091, 1449, 1477, 1967, 1677, 697, 1045, 737, 1117,
    1737, 667, 911, 1325, 473, 437, 1281, 1795, 1001, 261, 879, 51, 775, 1195, 801, 1635, 759, 165,
    1871, 1645, 1049, 245, 703, 1597, 553, 955, 209, 1779, 1849, 661, 865, 291, 841, 997, 1265,
    1965, 1625, 53, 1409, 893, 105, 1925, 1297, 589, 377, 1579, 929, 1053, 1655, 1829, 305, 1811,
    1895, 139, 575, 189, 343, 709, 1711, 1139, 1095, 277, 993, 1699, 55, 1435, 655, 1491, 1319,
    331, 1537, 515, 791, 507, 623, 1229, 1529, 1963, 1057, 355, 1545, 603, 1615, 1171, 743, 523,
    447, 1219, 1239, 1723, 465, 499, 57, 107, 1121, 989, 951, 229, 1521, 851, 167, 715, 1665, 1923,
    1687, 1157, 1553, 1869, 1415, 1749, 1185, 1763, 649, 1061, 561, 531, 409, 907, 319, 1469, 1961,
    59, 1455, 141, 1209, 491, 1249, 419, 1847, 1893, 399, 211, 985, 1099, 1793, 765, 1513, 1275,
    367, 1587, 263, 1365, 1313, 925, 247, 1371, 1359, 109, 1561, 1291, 191, 61, 1065, 1605, 721,
    781, 1735, 875, 1377, 1827, 1353, 539, 1777, 429, 1959, 1483, 1921, 643, 617, 389, 1809, 947,
    889, 981, 1441, 483, 1143, 293, 817, 749, 1383, 1675, 63, 1347, 169, 827, 1199, 1421, 583,
    1259, 1505, 861, 457, 1125, 143, 1069, 807, 1867, 2047, 2045, 279, 2043, 111, 307, 2041, 597,
    1569, 1891, 2039, 1957, 1103, 1389, 231, 2037, 65, 1341, 727, 837, 977, 2035, 569, 1643, 1633,
    547, 439, 1307, 2033, 1709, 345, 1845, 1919, 637, 1175, 379, 2031, 333, 903, 213, 1697, 797,
    1161, 475, 1073, 2029, 921, 1653, 193, 67, 1623, 1595, 943, 1395, 1721, 2027, 1761, 1955, 1335,
    357, 113, 1747, 1497, 1461, 1791, 771, 2025, 1285, 145, 973, 249, 171, 1825, 611, 265, 1189,
    847, 1427, 2023, 1269, 321, 1475, 1577, 69, 1233, 755, 1223, 1685, 1889, 733, 1865, 2021, 1807,
    1107, 1447, 1077, 1663, 1917, 1129, 1147, 1775, 1613, 1401, 555, 1953, 2019, 631, 1243, 1329,
    787, 871, 885, 449, 1213, 681, 1733, 687, 115, 71, 1301, 2017, 675, 969, 411, 369, 467, 295,
    693, 1535, 509, 233, 517, 401, 1843, 1543, 939, 2015, 669, 1527, 421, 591, 147, 281, 501, 577,
    195, 215, 699, 1489, 525, 1081, 917, 1951, 2013, 73, 1253, 1551, 173, 857, 309, 1407, 899, 663,
    1915, 1519, 1203, 391, 1323, 1887, 739, 1673, 2011, 1585, 493, 1433, 117, 705, 1603, 1111, 965,
    431, 1165, 1863, 533, 1823, 605, 823, 1179, 625, 813, 2009, 75, 1279, 1789, 1559, 251, 657,
    563, 761, 1707, 1759, 1949, 777, 347, 335, 1133, 1511, 267, 833, 1085, 2007, 1467, 1745, 1805,
    711, 149, 1695, 803, 1719, 485, 1295, 1453, 935, 459, 1151, 381, 1641, 1413, 1263, 77, 1913,
    2005, 1631, 541, 119, 1317, 1841, 1773, 359, 651, 961, 323, 1193, 197, 175, 1651, 441, 235,
    1567, 1885, 1481, 1947, 881, 2003, 217, 843, 1023, 1027, 745, 1019, 913, 717, 1031, 1621, 1503,
    867, 1015, 1115, 79, 1683, 793, 1035, 1089, 1731, 297, 1861, 2001, 1011, 1593, 619, 1439, 477,
    585, 283, 1039, 1363, 1369, 1227, 895, 1661, 151, 645, 1007, 1357, 121, 1237, 1375, 1821, 1911,
    549, 1999, 1043, 1945, 1419, 1217, 957, 599, 571, 81, 371, 1351, 1003, 1311, 931, 311, 1381,
    1137, 723, 1575, 1611, 767, 253, 1047, 1787, 1169, 1997, 1273, 853, 1247, 413, 1289, 1883, 177,
    403, 999, 1803, 1345, 451, 1495, 1093, 1839, 269, 199, 1387, 1183, 1757, 1207, 1051, 783, 83,
    423, 1995, 639, 1155, 1943, 123, 751, 1459, 1671, 469, 1119, 995, 393, 219, 1743, 237, 153,
    1909, 1473, 1859, 1705, 1339, 337, 909, 953, 1771, 1055, 349, 1993, 613, 1393, 557, 729, 1717,
    511, 1533, 1257, 1541, 1425, 819, 519, 85, 991, 1693, 503, 1445, 433, 877, 1305, 1525, 1601,
    829, 809, 325, 1583, 1549, 1991, 1941, 927, 1059, 1097, 1819, 527, 1197, 1881, 1333, 383, 125,
    361, 891, 495, 179, 633, 299, 863, 285, 1399, 987, 1487, 1517, 1639, 1141, 1729, 579, 87, 1989,
    593, 1907, 839, 1557, 799, 1629, 201, 155, 1649, 1837, 1063, 949, 255, 1283, 535, 773, 1681,
    461, 1785, 683, 735, 1123, 1801, 677, 689, 1939, 487, 757, 1857, 1987, 983, 443, 1327, 1267,
    313, 1173, 671, 221, 695, 1509, 271, 1619, 89, 565, 127, 1405, 1431, 1659, 239, 1101, 1159,
    1067, 607, 1565, 905, 1755, 1231, 1299, 665, 373, 1985, 701, 1879, 1221, 849, 627, 1465, 789,
    543, 1187, 1591, 923, 1905, 979, 1241, 181,
  ];
  println!("Step 2");

  // Step 2
  let copy = n.clone();
  let copy = (copy.clone() & BigUint::from(4294967295 as u64)) + (copy >> 32);
  let copy = (copy.clone() & BigUint::from(65535 as u16)) + (copy >> 16);
  let copy = (copy.clone() & BigUint::from(255 as u8))
    + ((copy.clone() >> 8) & BigUint::from(255 as u8))
    + (copy >> 16);
  // println!("{}", n.to_u64().unwrap());
  if bad255[copy.to_u64().unwrap() as usize] {
    return false;
  }

  println!("Step 3");

  let mut x = n.clone();
  if x.clone() & BigUint::from(4294967295 as u64) == zero {
    x >>= 32;
  }
  if x.clone() & BigUint::from(65535 as u64) == zero {
    x >>= 16;
  }
  if x.clone() & BigUint::from(255 as u16) == zero {
    x >>= 8;
  }
  if x.clone() & BigUint::from(15 as u8) == zero {
    x >>= 4;
  }
  if x.clone() & BigUint::from(3 as u8) == zero {
    x >>= 2;
  }
  if x.clone() & BigUint::from(7 as u8) != one {
    return false;
  }

  println!("Step 4");

  // let mut r: i64 = start[((n >> 3) & BigUint::from(1023 as u64)).to_u64().unwrap() as usize];
  // let mut t: BigInt;
  // let mut z: BigInt;
  // let zero_i = BigInt::from(0 as i8);
  // while {
  //   z = BigInt::from(x.clone()) - BigInt::from(r * r);
  //   if z == zero_i {
  //     return true;
  //   }
  //   t = z.clone() & -z.clone();
  //   r += ((z & t.clone()) >> 1).to_i64().unwrap();
  //   if r > (t.clone() >> 1).to_i64().unwrap() {
  //     r = t.to_i64().unwrap() - r;
  //   }
  //   t <= BigInt::from(1 << 33)
  // } {}
  // println!("All else fails");

  //0xC840C04048404040
  // let inbase16 = &[12, 8, 4, 0, 12, 0, 4, 0, 4, 8, 4, 0, 4, 0, 4, 0];
  // let good_mask = BigUint::from_radix_be(inbase16, 16).unwrap();
  // if good_mask << n >= zero {
  //   return false;
  // }
  true
}

#[allow(dead_code)]
fn passes_lucas(_n: &BigUint) -> bool {
  false
}

#[cfg(test)]
mod tests {
  use super::*;

  fn isqrt(n: usize) -> usize {
    n == 0 && return n;
    let mut s = (n as f64).sqrt() as usize;
    s = (s + n / s) >> 1;
    if s * s > n {
      s - 1
    } else {
      s
    }
  }

  fn perfect_sqrt(n: usize) -> isize {
    match n & 0xf {
      0 | 1 | 4 | 9 => {
        let t = isqrt(n);
        if t * t == n {
          t as isize
        } else {
          -1
        }
      }
      _ => -1,
    }
  }
  #[test]
  fn test_is_square() {
    let mut squares = vec![];
    for i in 0..200 {
      if perfect_sqrt(i) >= 0 {
        squares.push(true);
      } else {
        squares.push(false);
      }
    }
    for i in 0..squares.len() {
      println!("{}", i);
      let val = BigUint::from(i as u8);
      assert!(is_square(&val) == squares[i]);
    }
  }
  #[test]
  fn test_small_prime_factor() {
    let n_prime = BigUint::from(233u64);
    let n_composite = BigUint::from(50_621u64);
    let n_composite_large = BigUint::from(104_927u64);

    assert!(n_composite == BigUint::from(223u64) * BigUint::from(227u64));
    assert!(n_composite_large == BigUint::from(317u64) * BigUint::from(331u64));

    assert!(!has_small_prime_factor(&n_prime));
    assert!(has_small_prime_factor(&n_composite));
    assert!(!has_small_prime_factor(&n_composite_large));
  }

  #[test]
  fn test_miller_rabin() {
    assert!(passes_miller_rabin_base_2(&BigUint::from(13u64)));
    assert!(!passes_miller_rabin_base_2(&BigUint::from(65u64)));
  }
}
