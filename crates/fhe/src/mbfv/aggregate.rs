use crate::errors::Result;

/// Aggregate shares in an MPC protocol
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
        A::from_shares(iter.into_iter().collect::<Result<Vec<_>>>()?)
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
