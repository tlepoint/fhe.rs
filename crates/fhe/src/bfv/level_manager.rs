use std::sync::Arc;

use crate::{Error, Result};

use super::{context_chain::ContextLevel, Ciphertext};

/// Utility struct to manage ciphertext levels automatically
pub struct LevelManager {
    context_chain: Arc<ContextLevel>,
}

impl LevelManager {
    /// Create a new level manager from the head of a context chain
    pub fn new(context_chain: Arc<ContextLevel>) -> Self {
        Self { context_chain }
    }

    /// Determine the best level for addition
    pub fn optimal_level_for_add(&self, a: &Ciphertext, b: &Ciphertext) -> usize {
        a.level.max(b.level)
    }

    /// Determine the best level for multiplication
    pub fn optimal_level_for_mul(&self, a: &Ciphertext, b: &Ciphertext) -> Result<usize> {
        let next_level = self.optimal_level_for_add(a, b) + 1;
        if next_level <= self.context_chain.max_level() {
            Ok(next_level)
        } else {
            Err(Error::DefaultError(format!(
                "Invalid level: {next_level}, max level: {}",
                self.context_chain.max_level()
            )))
        }
    }

    /// Align two ciphertexts to the same level
    pub fn align_levels(&self, a: &mut Ciphertext, b: &mut Ciphertext) -> Result<()> {
        let target = self.optimal_level_for_add(a, b);
        a.switch_to_level(target)?;
        b.switch_to_level(target)?;
        Ok(())
    }

    /// Align a batch of ciphertexts
    pub fn align_batch(&self, cts: &mut [Ciphertext]) -> Result<usize> {
        if cts.is_empty() {
            return Ok(0);
        }
        let target = cts.iter().map(|c| c.level).max().unwrap();
        for ct in cts {
            ct.switch_to_level(target)?;
        }
        Ok(target)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::{BfvParametersBuilder, Ciphertext};

    #[test]
    fn align_two_ciphertexts() {
        let params = Arc::new(
            BfvParametersBuilder::new()
                .set_degree(16)
                .set_plaintext_modulus(1153)
                .set_moduli_sizes(&[50, 50])
                .build()
                .unwrap(),
        );

        let chain = params.context_chain();
        let manager = LevelManager::new(chain);

        let mut ct1 = Ciphertext::zero(&params);
        let mut ct2 = Ciphertext::zero(&params);
        ct2.switch_down().unwrap();

        manager.align_levels(&mut ct1, &mut ct2).unwrap();
        assert_eq!(ct1.level, ct2.level);
    }

    #[test]
    fn align_batch_ciphertexts() {
        let params = Arc::new(
            BfvParametersBuilder::new()
                .set_degree(16)
                .set_plaintext_modulus(1153)
                .set_moduli_sizes(&[50, 50])
                .build()
                .unwrap(),
        );

        let chain = params.context_chain();
        let manager = LevelManager::new(chain);

        let mut cts = vec![Ciphertext::zero(&params), Ciphertext::zero(&params)];
        cts[1].switch_down().unwrap();
        let level = manager.align_batch(&mut cts).unwrap();
        assert_eq!(level, cts[0].level);
        assert_eq!(cts[0].level, cts[1].level);
    }
}
