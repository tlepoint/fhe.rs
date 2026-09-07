use crate::errors::Result;

/// Aggregate shares in an MPC protocol.
///
/// The `Aggregate<Result<S>>` adapter forwards shares as they arrive without
/// retaining a vector. It stops on the first input error or aggregation error;
/// an observed input error takes precedence over an incomplete aggregate.
pub trait Aggregate<S>: Sized {
    /// Aggregate shares in an MPC protocol.
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = S>;
}

#[diagnostic::do_not_recommend]
impl<S, A> Aggregate<Result<S>> for A
where
    A: Aggregate<S>,
{
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = Result<S>>,
    {
        let mut input_error = None;
        let result = A::from_shares(
            iter.into_iter()
                .map_while(|share| match share {
                    Ok(share) => Some(share),
                    Err(error) => {
                        input_error = Some(error);
                        None
                    }
                })
                .fuse(),
        );
        match input_error {
            Some(error) => Err(error),
            None => result,
        }
    }
}

/// Perform aggregation directly on an iterator of shares.
///
/// This trait exists for convenience; the `aggregate` method is analogous to
/// [`Iterator::collect`], but the trait bound required is [`Aggregate`] rather
/// than [`FromIterator`].
pub trait AggregateIter {
    /// The type of share being aggregated.
    type Share;

    /// Aggregate shares in an MPC protocol.
    fn aggregate<A>(self) -> Result<A>
    where
        A: Aggregate<Self::Share>;
}

#[diagnostic::do_not_recommend]
impl<I: Iterator<Item = S>, S> AggregateIter for I {
    type Share = S;

    fn aggregate<A>(self) -> Result<A>
    where
        A: Aggregate<Self::Share>,
    {
        Aggregate::from_shares(self)
    }
}

#[cfg(test)]
mod tests {
    use super::{Aggregate, AggregateIter};
    use crate::errors::Result;

    #[derive(Debug, PartialEq, Eq)]
    struct Sum(u64);

    impl Aggregate<u64> for Sum {
        fn from_shares<T>(iter: T) -> Result<Self>
        where
            T: IntoIterator<Item = u64>,
        {
            Ok(Sum(iter.into_iter().sum()))
        }
    }

    #[test]
    fn aggregate_iter_collects_shares() -> Result<()> {
        let sum = vec![1u64, 2, 3].into_iter().aggregate::<Sum>()?;
        assert_eq!(sum, Sum(6));
        Ok(())
    }

    #[test]
    fn aggregate_result_flattens_shares() -> Result<()> {
        let sum = <Sum as Aggregate<Result<u64>>>::from_shares(vec![Ok(1u64), Ok(2), Ok(3)])?;
        assert_eq!(sum, Sum(6));
        Ok(())
    }
}

#[cfg(test)]
mod streaming_tests {
    use super::*;
    use std::cell::Cell;

    struct Count(usize);
    impl Aggregate<u8> for Count {
        fn from_shares<T: IntoIterator<Item = u8>>(iter: T) -> Result<Self> {
            Ok(Self(iter.into_iter().count()))
        }
    }

    struct PollAgain;
    impl Aggregate<u8> for PollAgain {
        fn from_shares<T: IntoIterator<Item = u8>>(iter: T) -> Result<Self> {
            let mut iter = iter.into_iter();
            assert_eq!(iter.next(), Some(1));
            assert_eq!(iter.next(), None);
            assert_eq!(iter.next(), None);
            Ok(Self)
        }
    }

    #[test]
    fn input_error_permanently_ends_the_adapter() {
        let visits = Cell::new(0);
        let result = PollAgain::from_shares(
            [
                Ok(1),
                Err(crate::MultipartyError::IncompatibleShares.into()),
                Ok(2),
            ]
            .into_iter()
            .inspect(|_| visits.set(visits.get() + 1)),
        );
        assert!(result.is_err());
        assert_eq!(visits.get(), 2);
    }

    #[test]
    fn result_adapter_stops_at_input_errors_without_collecting() {
        let visits = Cell::new(0);
        let inputs = [
            Ok(1),
            Err(crate::MultipartyError::IncompatibleShares.into()),
            Ok(2),
        ];
        let result =
            Count::from_shares(inputs.into_iter().inspect(|_| visits.set(visits.get() + 1)));
        assert!(matches!(
            result,
            Err(crate::Error::Multiparty(
                crate::MultipartyError::IncompatibleShares
            ))
        ));
        assert_eq!(visits.get(), 2);
        assert_eq!(Count::from_shares([Ok(1), Ok(2)]).unwrap().0, 2);
    }
}
