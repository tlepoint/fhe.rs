use crate::errors::Result;

/// Aggregate shares in an MPC protocol
pub trait Aggregate<S>: Sized {
    /// Aggregate shares in an MPC protocol.
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = S>;
}

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

impl<I: Iterator<Item = S>, S> AggregateIter for I {
    type Share = S;

    fn aggregate<A>(self) -> Result<A>
    where
        A: Aggregate<Self::Share>,
    {
        Aggregate::from_shares(self)
    }
}
